package main

import "fmt"
import "math"
import "math/rand"
import "flag"
import "strconv"
import "runtime"
import "sync"

// import "io/ioutil"
import "os"
import "bufio"

/** Vector Operations **/
type vector struct {
    x, y, z float64
}

func vector_add(vecA vector, vecB vector) vector {
    return vector{vecA.x + vecB.x, vecA.y + vecB.y, vecA.z + vecB.z}
}

func vector_sub(vecA vector, vecB vector) vector {
    return vector{vecA.x - vecB.x, vecA.y - vecB.y, vecA.z - vecB.z}
}

func vector_times(vec vector, c float64) vector {
    return vector{vec.x * c, vec.y * c, vec.z * c}
}

func vector_mul(vecA vector, vecB vector) vector {
    return vector{vecA.x * vecB.x, vecA.y * vecB.y, vecA.z * vecB.z}
}

func vector_dot(vecA vector, vecB vector) float64 { //cross
    return vecA.x*vecB.x + vecA.y*vecB.y + vecA.z*vecB.z
}

func vector_norm(vecA vector) vector {
    return vector_times(vecA, (1 / math.Sqrt(math.Pow(vecA.x, 2)+math.Pow(vecA.y, 2)+math.Pow(vecA.z, 2))))
}

func vector_cross(vecA vector, vecB vector) vector {
    //Vec(y*b.z-z*b.y,       z*b.x-x*b.z,          x*b.y-y*b.x)
    return vector{vecA.y*vecB.z - vecA.z*vecB.y, vecA.z*vecB.x - vecB.z, vecA.x*vecB.y - vecA.y*vecB.x}
}

/** END OF Vector operations **/
type ray struct {
    o vector
    d vector
}
type Refl_t int

const (
    DIFF = 1 << iota
    SPEC
    REFR
)

/** Sphere operations **/
type sphere struct {
    rad     float64 //radius
    p, e, c vector  // position, emission, color
    refl    int
}

func sphere_intersect(sphereA sphere, rayA ray) float64 {
    op := vector_sub(sphereA.p, rayA.o)
    var eps float64 = 1e-4
    var b float64 = vector_dot(op, rayA.d)
    var det float64 = b*b - vector_dot(op, op) + sphereA.rad * sphereA.rad
    if det < 0.0 {
        return 0.0
    } else {
        det = math.Sqrt(det)
    }
    if (b - det) > eps {
        return b - det
    } else if (b + det) > eps {
        return b + det
    }

    return 0.0
}

/** END OF Sphere operations **/
var spheres = []sphere{ //Scene: radius, position, emission, color, material
    sphere{1e5, vector{ 1e5+1,40.8,81.6}, vector{0.0,0.0,0.0},vector{.75,.25,.25},  DIFF},//Left
    sphere{1e5, vector{-1e5+99,40.8,81.6},vector{0.0,0.0,0.0},vector{.25,.55,.5},   REFR},//Rght
    sphere{1e5, vector{50,40.8, 1e5},     vector{0.0,0.0,0.0},vector{.75,.75,.75},  DIFF},//Back
    sphere{1e5, vector{50,40.8,-1e5+170}, vector{0.0,0.0,0.0},vector{0.0, 0.0, 0.0},DIFF},//Frnt
    sphere{1e5, vector{50, 1e5, 81.6},    vector{0.0,0.0,0.0},vector{.75,.75,.75},  DIFF},//Botm
    sphere{1e5, vector{50,-1e5+81.6,81.6},vector{0.0,0.0,0.0},vector{.75,.75,.75},  DIFF},//Top
    sphere{16.5,vector{27,16.5,47},       vector{0.0,0.0,0.0},vector{.999,.999,.999},SPEC},//Mirr
    sphere{16.5,vector{73,16.5,78},       vector{0.0,0.0,0.0},vector{.999,.999,.999},REFR},//Glas
    sphere{600, vector{50,681.6-.27,81.6},vector{12,12,12},   vector{0.0,0.0,0.0},   DIFF},//Lite
}

func clamp(x float64) float64 {
    if x < 0 {
        return 0
    } else if x > 1 {
        return 1
    }
    return x
}

func toByte(x float64) byte {
    return byte(math.Pow(clamp(x), 1/2.2) * 255 + 0.5)
}

func intersect(r *ray, t *float64, id *int) bool {
    var d, inf float64
    inf = 1e50
    *t = inf
    for i, s := range spheres {
        d = sphere_intersect(s,*r)
        if d != 0.0 && d < *t {
            *t = d
            *id = i
        }
    }
    return *t < inf
}

func radiance(r *ray, depth int) vector {
    var t float64                // distance to intersection
    id := 0                      // id of intersected object
    if !intersect(r, &t, &id) { // if miss, return black
        return vector{0, 0, 0}
    }

    obj := spheres[id] // the hit object

    if depth > 10 {
        return vector{0, 0, 0}
    }
    x := vector_add(r.o, vector_times(r.d, t)) // ray intersection point
    n := vector_norm(vector_sub(x, obj.p))     // sphere normal

    var nl vector
    if vector_dot(n, r.d) < 0.0 {
        nl = n
    } else {
        nl = vector_times(n, -1.0)
    }

    f := obj.c // object color (BRDF modulator)
    // use maximum reflectivity amount of Russian roulette
    var p float64
    if f.x > f.y && f.x > f.z {
        p = f.x
    } else if f.y > f.z {
        p = f.y
    } else {
        p = f.z
    }
    if depth > 4 {
        if erand48() < p {
            f = vector_times(f, 1.0/p)
        } else {
            return obj.e
        }
    }

    // ideal diffuse reflection
    if obj.refl == DIFF {
        var r1, r2, r2s float64
        r1 = 2 * math.Pi * rand.Float64()
        r2 = rand.Float64()
        r2s = math.Sqrt(r2)
        var w, u, v, d vector
        w = nl
        if math.Abs(w.x) > 0.1 {
            u = vector{0.0, 1.0, 0.0}
        } else {
            u = vector{1.0, 0.0, 0.0}
        }
        v = vector_cross(w, u)
        d = vector_norm(vector_add(vector_add(vector_times(u, math.Cos(r1) * r2s), vector_times(v, math.Sin(r1) * r2s)), vector_times(w, math.Sqrt(1 - r2))))
        return vector_add(obj.e, vector_mul(f, radiance(&ray{x, d}, depth+1)))
    } else if obj.refl == SPEC { // Ideal SPECULAR reflection

        return vector_add(obj.e, vector_mul(f, radiance(&ray{x, vector_sub(r.d, vector_times(n, 2 * vector_dot(n, r.d)))}, depth+1)))
    }
    var reflRay *ray = &ray{x, vector_sub(r.d, vector_times(n, 2*vector_dot(n, r.d)))}
    into := vector_dot(n, nl) > 0                                                            // Ray from outside going in?
    var nnt, ddn, cos2t float64
    nc := 1.0
    nt := 1.5
    if into {
        nnt = nc / nt
    } else {
        nnt = nt / nc
    }
    ddn = vector_dot(r.d, nl)
    cos2t = 1 - nnt * nnt * (1 - ddn * ddn)
    if cos2t < 0 { // Total internal reflection
        return vector_add(obj.e, vector_mul(f, radiance(reflRay, depth+1)))
    }

    var tdir vector  = vector_times(r.d, nnt)
    if into {
        tdir = vector_norm(vector_sub(tdir, vector_times(n, ddn*nnt+math.Sqrt(cos2t))))
    } else {
        tdir = vector_norm(vector_sub(tdir, vector_times(n, -(ddn*nnt+math.Sqrt(cos2t)))))
    }
    var a, b, c, R0, Re, Tr, P, RP, TP float64
    a, b = nt-nc, nt+nc
    R0 = a*a/(b*b)
    if into {
        c = 1 + ddn
    } else {
        c = 1 - vector_dot(tdir, n)
    }

    Re = R0 + (1 - R0) *c*c*c*c*c
    Tr = 1 - Re
    P = .25 + .5 * Re
    RP = Re/P
    TP = Tr/(1-P)

    if depth > 1 {
        if rand.Float64() < P {
            return vector_add(obj.e, vector_mul(f, vector_times(radiance(reflRay, depth+1), RP)))
        } else {
            return vector_add(obj.e, vector_mul(f, vector_times(radiance(&ray{x, tdir}, depth+1), TP)))
        }
    }
    return vector_add(obj.e, vector_mul(f, vector_times(vector_add(radiance(reflRay, depth+1),
        radiance(&ray{x, tdir}, depth+1)), Tr)))
}

func check(e error) {
    if e != nil {
        panic(e)
    }
}

func erand48() float64 {
    return rand.Float64()
}

type job_input struct {
    x int
    y int
}

type job_result struct {
    radiance vector
    i int
}

func main() {
    var NUM_PROCS = runtime.NumCPU()
    runtime.GOMAXPROCS(NUM_PROCS)
    flag.Parse()
    w := 1024
    h := 768
    samps := 1

    if flag.NArg() > 0 {
        samps, _ = strconv.Atoi(flag.Arg(0))
        samps /= 4
    }

    // Cam Setup
    cam := &ray{vector{50, 52, 295.6}, vector_norm(vector{0, -0.042612, -1})} // cam pos, dir
    cx := vector{float64(w) * 0.5135 / float64(h), 0, 0}
    cy := vector_times(vector_cross(cx, cam.d), 0.5135)
    c := make([]vector, w*h)



    jobs := make(chan job_input)
    result_channel := make(chan job_result)

    var wg sync.WaitGroup
    wg.Add(NUM_PROCS)

    for i := 0; i < NUM_PROCS; i++ {
        go func() {
            var direction, r vector
            var r1, r2 float64
            var dx, dy float64
            for {
                j, more := <-jobs
                if more {
                    //fmt.Printf("job input (%d,%d) \n", j.x, j.y)
                    fmt.Printf("\rRendering (%d spp) %5.2f",samps*4,100.*float64(j.y)/(float64(h)-1.0));
                    for sy := 0.0; sy < 2.0; sy += 1.0 { // 2x2 subpixel rows
                        for sx := 0.0; sx < 2.0; sx += 1.0 { // 2x2 subpixel cols
                            r = vector{0, 0, 0}
                            for s := 0; s < samps; s++ {
                                r1, r2 = 2 * erand48(), 2 * erand48()
                                if r1 < 1 {
                                    dx = math.Sqrt(r1) - 1
                                } else {
                                    dx = 1 - math.Sqrt(2-r1)
                                }
                                if r2 < 1 {
                                    dy = math.Sqrt(r2) -1
                                } else {
                                    dy = 1 - math.Sqrt(2-r2)
                                }

                                direction = vector_add(vector_add(vector_times(cx, ((float64(sx)*.5+dx)/2+float64(j.x))/float64(w)-.5),
                                    vector_times(cy, ((float64(sy)+.5+dy)/2+float64(j.y))/float64(h)-.5)), cam.d)
                                r = vector_add(r, vector_times(radiance(&ray{vector_add(cam.o, vector_times(direction, 140.0)), vector_norm(direction)}, 0), 1.0/float64(samps)))
                            } // Camera rays are pushed ^^^^^ forward to start in interior
                            if r.x != 0 || r.y !=0 || r.z != 0 {
                                res_i := (h-j.y-1)*w+j.x
                                result_channel <- job_result{r, res_i}
                            }
                        }
                    }
                } else {
                    fmt.Println("received all jobs")
                    defer wg.Done()
                    return
                }
            }
        }()
    }

    go func() {
        for {
            job_res, more := <-result_channel
            if more {
                // fmt.Println(job_res)
                c[job_res.i] = vector_add(c[job_res.i], vector_times(job_res.radiance, 0.25))
            } else {
                fmt.Println("Closing result job result goroutine.")
                return
            }
        }
    }()

    for y := 0; y < h; y++ { // Loop over image rows
        for x := 0; x < w; x++ { // Loop cols
            jobs <- job_input{x, y}
        }
    }
    close(jobs)
    wg.Wait()


    f, err := os.Create("image.ppm")
    check(err)
    wfile := bufio.NewWriter(f)

    fmt.Fprintf(wfile, "P3\n%d %d\n%d\n", w, h, 255)
    for i := 0; i < w*h; i++ {
        fmt.Fprintf(wfile, "%d %d %d ", toByte(c[i].x), toByte(c[i].y), toByte(c[i].z))
    }

    wfile.Flush()
}

/*
    vec1 := vector{-1, -3, -5}
    vec2 := vector{2, 2, 2}
    fmt.Println(vector_add(vec1, vec2))
    fmt.Println(vector_sub(vec1, vec2))
    fmt.Println(vector_times(vec1, 4))
    fmt.Println(vector_mul(vec1, vec2))
    fmt.Println(vector_dot(vec1, vec2))
    fmt.Println(vector_norm(vec1))
    fmt.Println(vector_cross(vec1, vec2))

    f, err := os.Create("/tmp/dat2")
    check(err)
    w := bufio.NewWriter(f)

    fmt.Fprintln(w, "This is a test.")
    w.Flush()
*/
