// Package hermit3 contains a float32 type T and functions for 3D cubic hermit splines.
// See: http://en.wikipedia.org/wiki/Cubic_Hermite_spline
package hermit

import (
	"fmt"
	math64 "math"

	"github.com/ungerik/go3d/vec3"
)

// PointTangent contains a point and a tangent at that point.
// This is a helper sub-struct for T.
type PointTangent struct {
	Point   vec3.T
	Tangent vec3.T
}

// T holds the data to define a hermit spline.
type T struct {
	A PointTangent
	B PointTangent
}

// Parse parses T from a string. See also String()
func Parse(s string) (r T, err error) {
	_, err = fmt.Sscan(s,
		&r.A.Point[0], &r.A.Point[1], &r.A.Point[2],
		&r.A.Tangent[0], &r.A.Tangent[1], &r.A.Tangent[2],
		&r.B.Point[0], &r.B.Point[1], &r.B.Point[2],
		&r.B.Tangent[0], &r.B.Tangent[1], &r.B.Tangent[2],
	)
	return r, err
}

// String formats T as string. See also Parse().
func (herm *T) String() string {
	return fmt.Sprintf("%s %s %s %s",
		herm.A.Point.String(), herm.A.Tangent.String(),
		herm.B.Point.String(), herm.B.Tangent.String(),
	)
}

// Point returns a point on a hermit spline at t (0,1).
func (herm *T) Point(t float32) vec3.T {
	return Point(&herm.A.Point, &herm.A.Tangent, &herm.B.Point, &herm.B.Tangent, t)
}

// Tangent returns the derivative dPoint/dt of the hermit spline at t (0,1).
// The result is not normalized; its length is the curve speed at t.
// At t=0 it equals A.Tangent and at t=1 it equals B.Tangent, as defined by the
// cubic Hermite basis. If any of the four control vectors is not finite,
// Tangent returns a vector with all components set to NaN.
func (herm *T) Tangent(t float32) vec3.T {
	return Tangent(&herm.A.Point, &herm.A.Tangent, &herm.B.Point, &herm.B.Tangent, t)
}

// Length returns the arc length of the hermit spline from t=0 to t.
// The result is non-negative and non-decreasing in t. Values of t outside
// (0,1) are clamped to the spline domain [0,1], and a NaN t yields NaN.
// If any of the four control vectors is not finite, Length returns NaN.
// A curve whose derivative is identically zero has length zero.
func (herm *T) Length(t float32) float32 {
	return Length(&herm.A.Point, &herm.A.Tangent, &herm.B.Point, &herm.B.Tangent, t)
}

// Point returns a point on a hermit spline at t (0,1).
func Point(pointA, tangentA, pointB, tangentB *vec3.T, t float32) vec3.T {
	t2 := t * t
	t3 := t2 * t

	f := 2*t3 - 3*t2 + 1
	result := pointA.Scaled(f)

	f = t3 - 2*t2 + t
	tAf := tangentA.Scaled(f)
	result.Add(&tAf)

	f = t3 - t2
	tBf := tangentB.Scaled(f)
	result.Add(&tBf)

	f = -2*t3 + 3*t2
	pAf := pointB.Scaled(f)
	result.Add(&pAf)

	return result
}

// Tangent returns the derivative dPoint/dt of the cubic Hermite spline at t.
// The result is not normalized: its length is the speed of the curve at t.
// By construction Tangent(0) equals tangentA and Tangent(1) equals tangentB.
// Values of t outside [0,1] extrapolate the cubic polynomial, matching Point.
// If any control vector is not finite, or t is NaN, the result is all NaN.
func Tangent(pointA, tangentA, pointB, tangentB *vec3.T, t float32) vec3.T {
	if !finiteControls3(pointA, tangentA, pointB, tangentB) || !isFinite64(float64(t)) {
		nan := float32(math64.NaN())
		return vec3.T{nan, nan, nan}
	}

	tt := float64(t)
	t2 := tt * tt

	// Derivatives of the cubic Hermite basis functions:
	// h00 = 2t^3 - 3t^2 + 1, h10 = t^3 - 2t^2 + t,
	// h01 = -2t^3 + 3t^2,    h11 = t^3 - t^2.
	d00 := 6*t2 - 6*tt
	d10 := 3*t2 - 4*tt + 1
	d01 := -6*t2 + 6*tt
	d11 := 3*t2 - 2*tt

	var result vec3.T
	for i := 0; i < 3; i++ {
		result[i] = float32(float64(pointA[i])*d00 +
			float64(tangentA[i])*d10 +
			float64(tangentB[i])*d11 +
			float64(pointB[i])*d01)
	}
	return result
}

// Length returns the arc length of the cubic Hermite spline from t=0 to t.
// The result is non-negative and non-decreasing in t: values of t <= 0 give 0,
// values of t > 1 are clamped to 1, and a NaN t gives NaN. Curves whose
// derivative vanishes everywhere have length 0.
//
// The integral of |dPoint/dt| is evaluated with adaptive Simpson quadrature
// in float64 precision, so the result is not limited by a fixed sampling
// resolution. If any control vector is not finite, the result is NaN.
func Length(pointA, tangentA, pointB, tangentB *vec3.T, t float32) float32 {
	if !finiteControls3(pointA, tangentA, pointB, tangentB) || !isFinite64(float64(t)) {
		return float32(math64.NaN())
	}

	tt := float64(t)
	if tt <= 0 {
		return 0
	}
	if tt > 1 {
		tt = 1
	}

	speed := func(x float64) float64 {
		return speed3(pointA, tangentA, pointB, tangentB, x)
	}
	scale := 1 + curveScale3(pointA, tangentA, pointB, tangentB)
	return float32(adaptiveSimpson(speed, 0, tt, arcLengthTolerance*scale))
}

// speed3 returns |dPoint/dt| at t for a 3D cubic Hermite spline in float64.
func speed3(pointA, tangentA, pointB, tangentB *vec3.T, t float64) float64 {
	t2 := t * t
	d00 := 6*t2 - 6*t
	d10 := 3*t2 - 4*t + 1
	d01 := -6*t2 + 6*t
	d11 := 3*t2 - 2*t

	dx := float64(pointA[0])*d00 + float64(tangentA[0])*d10 +
		float64(tangentB[0])*d11 + float64(pointB[0])*d01
	dy := float64(pointA[1])*d00 + float64(tangentA[1])*d10 +
		float64(tangentB[1])*d11 + float64(pointB[1])*d01
	dz := float64(pointA[2])*d00 + float64(tangentA[2])*d10 +
		float64(tangentB[2])*d11 + float64(pointB[2])*d01

	return math64.Hypot(math64.Hypot(dx, dy), dz)
}

// curveScale3 returns a characteristic magnitude of the control data,
// used to scale the absolute quadrature tolerance.
func curveScale3(pointA, tangentA, pointB, tangentB *vec3.T) float64 {
	return math64.Hypot(math64.Hypot(float64(pointA[0]), float64(pointA[1])), float64(pointA[2])) +
		math64.Hypot(math64.Hypot(float64(tangentA[0]), float64(tangentA[1])), float64(tangentA[2])) +
		math64.Hypot(math64.Hypot(float64(tangentB[0]), float64(tangentB[1])), float64(tangentB[2])) +
		math64.Hypot(math64.Hypot(float64(pointB[0]), float64(pointB[1])), float64(pointB[2]))
}

// finiteControls3 reports whether every component of the control vectors
// is finite (neither NaN nor infinity).
func finiteControls3(pointA, tangentA, pointB, tangentB *vec3.T) bool {
	return finiteVec3(pointA) && finiteVec3(tangentA) &&
		finiteVec3(pointB) && finiteVec3(tangentB)
}

func finiteVec3(v *vec3.T) bool {
	for i := 0; i < 3; i++ {
		if !isFinite64(float64(v[i])) {
			return false
		}
	}
	return true
}

// isFinite64 reports whether f is neither NaN nor an infinity.
func isFinite64(f float64) bool {
	return !math64.IsNaN(f) && !math64.IsInf(f, 0)
}

// arcLengthTolerance is the relative accuracy target of the adaptive
// quadrature. It is far below float32 resolution so that rounding to the
// float32 result cannot make Length decrease between increasing t values.
const arcLengthTolerance = 1e-11

// maxSimpsonDepth bounds the adaptive subdivision, guaranteeing termination
// even for curves with cusps where the speed is not differentiable.
const maxSimpsonDepth = 40

// adaptiveSimpson integrates f over [a,b] to absolute accuracy eps.
func adaptiveSimpson(f func(float64) float64, a, b, eps float64) float64 {
	if b <= a {
		return 0
	}
	m := 0.5 * (a + b)
	fa, fm, fb := f(a), f(m), f(b)
	return simpsonStep(f, a, b, fa, fm, fb, eps, maxSimpsonDepth)
}

func simpsonStep(f func(float64) float64, a, b, fa, fm, fb, eps float64, depth int) float64 {
	m := 0.5 * (a + b)
	lm := 0.5 * (a + m)
	rm := 0.5 * (m + b)
	flm := f(lm)
	frm := f(rm)

	left := (m - a) / 6 * (fa + 4*flm + fm)
	right := (b - m) / 6 * (fm + 4*frm + fb)
	whole := (b - a) / 6 * (fa + 4*fm + fb)
	delta := left + right - whole

	if depth <= 0 || math64.Abs(delta) <= 15*eps {
		return left + right + delta/15
	}
	return simpsonStep(f, a, m, fa, flm, fm, eps/2, depth-1) +
		simpsonStep(f, m, b, fm, frm, fb, eps/2, depth-1)
}
