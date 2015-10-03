package main

import "fmt"
import "math"
import "math/rand"

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

func vector_mod(vecA vector, vecB vector) vector {
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
	eps := 10 * math.Exp(-4)
	b := vector_dot(op, rayA.d) // returns float
	det := b*b - vector_dot(op, op) + sphereA.rad*sphereA.rad

	if det < 0 {
		return 0
	} else {
		det = math.Sqrt(det)
		t := b - det
		if t > eps {
			return t
		} else {
			t := b + det
			if t > eps {
				return t
			} else {
				return 0
			}
		}
	}
}

/** END OF Sphere operations **/
var spheres = []sphere{ //Scene: radius, position, emission, color, material
	sphere{math.Exp(5), vector{math.Exp(5) + 1, 40.8, 81.6}, vector{0, 0, 0}, vector{.75, .25, .25}, DIFF},         //Left
	sphere{1 * math.Exp(5), vector{-1*math.Exp(5) + 99, 40.8, 81.6}, vector{0, 0, 0}, vector{.25, .25, .75}, DIFF}, //Rght
	sphere{1 * math.Exp(5), vector{50, 40.8, 1 * math.Exp(5)}, vector{0, 0, 0}, vector{.75, .75, .75}, DIFF},       //Back
	sphere{1 * math.Exp(5), vector{50, 40.8, -1*math.Exp(5) + 170}, vector{0, 0, 0}, vector{0, 0, 0}, DIFF},        //Frnt
	sphere{1 * math.Exp(5), vector{50, 1 * math.Exp(5), 81.6}, vector{0, 0, 0}, vector{.75, .75, .75}, DIFF},       //Botm
	sphere{1 * math.Exp(5), vector{50, -1*math.Exp(5) + 81.6, 81.6}, vector{0, 0, 0}, vector{.75, .75, .75}, DIFF}, //Top
	sphere{16.5, vector{27, 16.5, 47}, vector{0, 0, 0}, vector_times(vector{1, 1, 1}, 0.999), SPEC},                //Mirr
	sphere{16.5, vector{73, 16.5, 78}, vector{0, 0, 0}, vector_times(vector{1, 1, 1}, 0.999), REFR},                //Glass
	sphere{600, vector{50, 681.6 - .27, 81.6}, vector{12, 12, 12}, vector{0, 0, 0}, DIFF},                          //Lite
}

func clamp(x float64) float64 {
	if x < 0 {
		return 0
	}
	if x > 1 {
		return 1
	}
	return x
}

func toInt(x float64) int {
	return int(math.Pow(clamp(x), 1/2.2)*255 + .5)
}

func intersect(r *ray, t *float64, id *int) bool {
	n := len(spheres)
	d := 0.0
	inf := 1 * math.Exp(20)
	*t = inf
	for i := n - 1; i > 0; i-- {
		d = sphere_intersect(spheres[i], *r)
		if (d > 0) && (d < *t) {
			*t = d
			*id = i
		}
	}
	return *t < inf
}

func radiance(r ray, depth int, Xi [3]int) vector {
	var t float64                // distance to intersection
	id := 0                      // id of intersected object
	if !intersect(&r, &t, &id) { // if miss, return black
		return vector{0, 0, 0}
	}

	obj := spheres[id] // the hit object

	if depth > 10 {
		return vector{0, 0, 0}
	}
	x := vector_add(r.o, vector_times(r.d, t)) // ray intersection point
	n := vector_norm(vector_sub(x, obj.p))     // sphere normal
	nl := vector_times(n, -1)
	if vector_dot(n, (r.d)) < 0 { // properly oriented surface normal
		nl = n
	}

	f := obj.c // object color (BRDF modulator)
	// use maximum reflectivity amount of Russian roulette
	var p float64
	if f.x > f.y && f.x > f.z { // max refl
		p = f.x
	} else {
		if f.y > f.z {
			p = f.y
		} else {
			p = f.z
		}
	}

	depth++
	if depth > 5 || p == 0 {
		if erand48() < p {
			f = vector_times(f, (1 / p))
		} else {
			return obj.e //R.R.
		}
	}
	// ideal diffuse reflection
	if obj.refl == DIFF {
		r1 := 2 * math.Pi * erand48() // angle around
		r2 := erand48()
		r2s := math.Sqrt(r2) // distance from center
		w := nl
		// u is perpendicular to w:
		var u vector
		if math.Abs(w.x) > .1 {
			u = vector{0, 1, 0}
		} else {
			u = vector_mod(vector{1, 0, 0}, w)
		}
		u = vector_norm(u)
		v := vector_mod(w, u) // v is perpendicular to u and w
		d := vector_times(u, math.Cos(r1)*r2s)
		d = vector_add(d, vector_times(v, math.Sin(r1)*r2s))
		d = vector_add(d, vector_norm(vector_times(w, math.Sqrt(1-r2)))) // d is random reflection ray
		return vector_add(obj.e, vector_mul(f, (radiance(ray{x, d}, depth, Xi))))
	} else if obj.refl == SPEC { // Ideal SPECULAR reflection
		_ray := ray{x, vector_times(vector_times(vector_sub(r.d, n), 2), vector_dot(n, r.d))}
		return vector_add(obj.e, vector_mul(f, (radiance(_ray, depth, Xi)))) //fix this..
	}
	// Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION
	reflRay := ray{x, vector_sub(r.d, vector_times(vector_times(n, 2), vector_dot(n, r.d)))} // Ideal dielectric REFRACTION
	into := vector_dot(n, nl) > 0                                                            // Ray from outside going in?
	nc := 1.0
	nt := 1.5
	nnt := nt / nc
	if into {
		nnt = nc / nt
	}
	ddn := vector_dot(r.d, nl)
	cos2t := 1 - nnt*nnt*(1-ddn*ddn)
	if cos2t < 0 { // Total internal reflection
		return vector_add(obj.e, vector_mul(f, radiance(reflRay, depth, Xi)))
	}

	into_mul := -1.0
	if into {
		into_mul = 1.0
	}
	tdir := vector_norm(vector_sub(vector_times(r.d, nnt), vector_times(vector_times(n, into_mul), (ddn*nnt+math.Sqrt(cos2t)))))
	a := nt - nc
	b := nt + nc
	R0 := a * a / (b * b)

	c := 1 - (vector_dot(tdir, n))
	if into {
		c = 1 + ddn
	}
	//double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
	Re := R0 + (1-R0)*math.Pow(c, 5)
	Tr := 1 - Re
	P := 0.25 + .5*Re
	RP := Re / P
	TP := Tr / (1 - P)
	/* return obj.e + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette
	   radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
	*/
	if depth > 2 {
		if erand48() < P { // Russian roulette
			return vector_add(obj.e, vector_times(radiance(reflRay, depth, Xi), RP))
		} else {
			return vector_add(obj.e, vector_times(radiance(ray{x, tdir}, depth, Xi), TP))
		}
	} else {
		return vector_add(obj.e, vector_times(vector_add(vector_times(radiance(reflRay, depth, Xi), Re), radiance(ray{x, tdir}, depth, Xi)), Tr))
	}
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func erand48() float64 {
	return rand.Float64()
}

func main() {
	w := 1024
	h := 768
	samps := 100

	cam := ray{vector{50, 52, 295.6}, vector_norm(vector{0, -0.042612, -1})} // cam pos, dir
	cx := vector{float64(w) * 0.5135 / float64(h), 0, 0}
	cy := vector_times(vector_norm(vector_mod(cx, cam.d)), 0.5135)
	c := make([]vector, w*h)
	for y := 0; y < h; y++ { // Loop over image rows
		fmt.Printf("\rRendering (%d spp) %5.2f%%", samps*4, 100.*y/(h-1))
		Xi := [3]int{0, 0, y * y * y}
		for x := 0; x < w; x++ { // Loop cols
			i := (h-y-1)*w + x
			for sy := 0.0; sy < 2.0; sy += 1.0 { // 2x2 subpixel rows
				for sx := 0.0; sx < 2.0; sx += 1.0 { // 2x2 subpixel cols
					r := vector{0, 0, 0}
					for s := 0; s < samps; s++ {
						r1 := 2 * erand48()
						dx := 1 - math.Sqrt(2-r1)
						if r1 < 1 {
							dx = math.Sqrt(r1) - 1
						}

						r2 := 2 * erand48()
						dy := 1 - math.Sqrt(2-r2)
						if r2 < 1 {
							dy = math.Sqrt(r2) - 1
						}

						d := vector_add(cam.d, vector_add(vector_times(cx, (((sx+.5+dx)/2.0+float64(x))/float64(w)-.5)), vector_times(cy, (((float64(sy)+.5+dy)/2.0+float64(y))/float64(h)-.5))))

						r = vector_add(r, vector_times(radiance(ray{vector_add(cam.o, vector_times(d, 140)), vector_norm(d)}, 0, Xi), (1./float64(samps))))
					} // Camera rays are pushed ^^^^^ forward to start in interior
					c[i] = vector_times(vector_add(c[i], vector{clamp(r.x), clamp(r.y), clamp(r.z)}), .25)
				}
			}
		}
	}

	f, err := os.Create("image.ppm")
	check(err)
	wfile := bufio.NewWriter(f)

	fmt.Fprintf(wfile, "P3\n%d %d\n%d\n", w, h, 255)
	for i := 0; i < w*h; i++ {
		c[i] = vector_times(c[i], 255)
		fmt.Fprintf(wfile, "%d %d %d ", int(c[i].x), int(c[i].y), int(c[i].z))
		// fmt.Printf("%f %f %f \n", c[i].x, c[i].y, c[i].z)
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
	fmt.Println(vector_mod(vec1, vec2))

	f, err := os.Create("/tmp/dat2")
	check(err)
	w := bufio.NewWriter(f)

	fmt.Fprintln(w, "This is a test.")
	w.Flush()
*/
