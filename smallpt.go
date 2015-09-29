package main

import "fmt"
import "math"

/** Vector Operations **/
type vector struct {
    x, y, z float64
}
func vector_add(vecA vector, vecB vector) vector {
    return vector{vecA.x+vecB.x, vecA.y+vecB.y, vecA.z+vecB.z}
}
func vector_sub(vecA vector, vecB vector) vector {
    return vector{vecA.x-vecB.x, vecA.y-vecB.y, vecA.z-vecB.z}
}
func vector_times(vec vector, c float64) vector {
    return vector{vec.x*c, vec.y*c, vec.z*c}
}
func vector_mul(vecA vector, vecB vector) vector {
    return vector{vecA.x*vecB.x, vecA.y*vecB.y, vecA.z*vecB.z}
}
func vector_dot(vecA vector, vecB vector) float64 { //cross
    return vecA.x*vecB.x + vecA.y*vecB.y + vecA.z*vecB.z
}
func vector_norm(vecA vector) float64 {
    return 1/math.Sqrt(math.Pow(vecA.x,2)+math.Pow(vecA.y,2)+math.Pow(vecA.z,2))
}
func vector_mod(vecA vector, vecB vector) vector {
    //Vec(y*b.z-z*b.y,       z*b.x-x*b.z,          x*b.y-y*b.x)
    return vector{vecA.y*vecB.z-vecA.z*vecB.y, vecA.z*vecB.x-vecB.z, vecA.x*vecB.y-vecA.y*vecB.x}
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
    rad float64 //radius
    p, e, c vector // position, emission, color
    refl int
}

func sphere_intersect(sphereA sphere, rayA ray) float64 {
    op := vector_sub(sphereA.p, rayA.o)
    eps := 10*math.Exp(-4)
    b := vector_dot(op, rayA.d) // returns float
    det := b*b - vector_dot(op, op) + sphereA.rad*sphereA.rad

    if det < 0 {
        return 0
    } else {
        det = math.Sqrt(det);      
        t := b - det
        if t > eps {
            return t
        } else {
            t := b+det
            if t > eps {
                return t
            } else {
                return 0
            }
        }
    }
}
/** END OF Sphere operations **/
var spheres = []sphere{//Scene: radius, position, emission, color, material
  sphere{math.Exp(5), vector{math.Exp(5)+1,40.8,81.6}, vector{0,0,0}, vector{.75,.25,.25},     DIFF},//Left
  sphere{1*math.Exp(5), vector{-1*math.Exp(5)+99,40.8,81.6},vector{0,0,0},vector{.25,.25,.75}, DIFF},//Rght
  sphere{1*math.Exp(5), vector{50,40.8, 1*math.Exp(5)},     vector{0,0,0},vector{.75,.75,.75}, DIFF},//Back
  sphere{1*math.Exp(5), vector{50,40.8,-1*math.Exp(5)+170}, vector{0,0,0},vector{0,0,0},       DIFF},//Frnt
  sphere{1*math.Exp(5), vector{50, 1*math.Exp(5), 81.6},    vector{0,0,0},vector{.75,.75,.75}, DIFF},//Botm
  sphere{1*math.Exp(5), vector{50,-1*math.Exp(5)+81.6,81.6},vector{0,0,0},vector{.75,.75,.75}, DIFF},//Top
  sphere{16.5, vector{27,16.5,47}, vector{0,0,0}, vector_times(vector{1,1,1},0.999), SPEC},    //Mirr
  sphere{16.5,vector{73,16.5,78},  vector{0,0,0}, vector_times(vector{1,1,1},0.999), REFR},    //Glass
  sphere{600, vector{50,681.6-.27,81.6}, vector{12,12,12},  vector{0,0,0}, DIFF},     //Lite
}

func clamp(x float64 ) float64 { 
    if x < 0 {
        return 0
    } 
    if x > 1 {
        return 1
    } 
    return x
}

func toInt(x float64) int {
    return int(math.Pow(clamp(x),1/2.2)*255+.5);
}

func intersect(r *ray, t *float64, id *int) bool {
    n := len(spheres)
    d := 0.0
    inf := 1*math.Exp(20)
    *t = inf
    for i := n; i > 0; i-- {
        d = sphere_intersect(spheres[i], *r)
        if (d > 0) && (d < *t) {
            *t = d
            *id = i
        }
    }
    return *t < inf
}

func main() {
    // create spheres.
    /*

    */

    vec1 := vector{-1,-3,-5}
    vec2 := vector{2,2,2}
    fmt.Println(vector_add(vec1, vec2))
    fmt.Println(vector_sub(vec1, vec2))
    fmt.Println(vector_times(vec1, 4))
    fmt.Println(vector_mul(vec1, vec2))
    fmt.Println(vector_dot(vec1, vec2))
    fmt.Println(vector_norm(vec1))
    fmt.Println(vector_mod(vec1, vec2))
}
