package hermit

import (
	"errors"
	"math"

	"github.com/ungerik/go3d/vec3"
)

// ErrNonFiniteControl is returned by TryBoundingBox when a control point or
// tangent contains a NaN or infinite component.
var ErrNonFiniteControl = errors.New("hermit3: control point or tangent must not be NaN or Inf")

// BoundingBox returns the smallest axis aligned vec3.Box that contains every
// point of the Hermite segment for the whole parameter interval t in [0,1].
//
// The min/max of each component are found analytically from the endpoints and
// the (up to two) interior roots of the derivative of the cubic, so local
// extrema are not missed by sampling a fixed finite set of parameters.
//
// Extrema are evaluated in float64 and the resulting bounds are rounded
// outwards to float32, so the returned box conservatively contains the true
// curve while staying tight to roughly one float32 ULP.
//
// BoundingBox panics if any component of pointA, tangentA, pointB or tangentB
// is NaN or infinite. Use TryBoundingBox to receive the failure as an error.
//
// The receiver is a value so the query can be called directly on a composite
// literal such as T{...}.BoundingBox().
func (herm T) BoundingBox() vec3.Box {
	return BoundingBox(&herm.A.Point, &herm.A.Tangent, &herm.B.Point, &herm.B.Tangent)
}

// TryBoundingBox is like BoundingBox but returns ErrNonFiniteControl when a
// control value is non-finite instead of panicking.
func (herm T) TryBoundingBox() (vec3.Box, error) {
	return TryBoundingBox(&herm.A.Point, &herm.A.Tangent, &herm.B.Point, &herm.B.Tangent)
}

// Bounds is an alias for BoundingBox, matching the common Bounds naming.
func (herm T) Bounds() vec3.Box {
	return herm.BoundingBox()
}

// BoundingBox returns the smallest axis aligned vec3.Box containing the whole
// Hermite segment defined by pointA/tangentA and pointB/tangentB for t in [0,1].
//
// It panics on non-finite control values; see TryBoundingBox for an error
// returning variant.
func BoundingBox(pointA, tangentA, pointB, tangentB *vec3.T) vec3.Box {
	box, err := TryBoundingBox(pointA, tangentA, pointB, tangentB)
	if err != nil {
		panic(err)
	}
	return box
}

// TryBoundingBox is the error returning variant of BoundingBox.
func TryBoundingBox(pointA, tangentA, pointB, tangentB *vec3.T) (vec3.Box, error) {
	if err := validateControls(pointA, tangentA, pointB, tangentB); err != nil {
		return vec3.Box{}, err
	}

	var box vec3.Box
	for axis := 0; axis < 3; axis++ {
		lo, hi := hermiteAxisBounds64(
			float64(pointA[axis]), float64(tangentA[axis]),
			float64(pointB[axis]), float64(tangentB[axis]),
		)
		// Round the exact float64 extrema outwards to the nearest float32 so
		// that the float32 box does not dip inside the true curve. Exact
		// float32 values (such as the canonical 0 and 1) round unchanged.
		min32 := float32(lo)
		if float64(min32) > lo {
			min32 = math.Nextafter32(min32, float32(math.Inf(-1)))
		}
		max32 := float32(hi)
		if float64(max32) < hi {
			max32 = math.Nextafter32(max32, float32(math.Inf(1)))
		}
		box.Min[axis] = min32
		box.Max[axis] = max32
	}
	return box, nil
}

// Bounds is an alias for BoundingBox.
func Bounds(pointA, tangentA, pointB, tangentB *vec3.T) vec3.Box {
	return BoundingBox(pointA, tangentA, pointB, tangentB)
}

func validateControls(pointA, tangentA, pointB, tangentB *vec3.T) error {
	controls := []struct {
		name string
		v    *vec3.T
	}{
		{"pointA", pointA},
		{"tangentA", tangentA},
		{"pointB", pointB},
		{"tangentB", tangentB},
	}
	for _, c := range controls {
		if c.v == nil {
			return errors.New("hermit3: " + c.name + " must not be nil")
		}
		for axis := 0; axis < 3; axis++ {
			f := float64(c.v[axis])
			if math.IsNaN(f) || math.IsInf(f, 0) {
				return errors.New("hermit3: " + c.name + " contains a NaN or Inf component")
			}
		}
	}
	return nil
}

// hermiteAxisBounds64 returns the float64 minimum and maximum of one component
// of the cubic Hermite curve over t in [0,1]. Control values originate from
// float32 but are converted exactly, which keeps root finding and extrema
// evaluation stable for cancellation and extremely small scales.
//
// The extrema are the endpoint values and the values at any interior root of
// the derivative, solved analytically with a numerically stable quadratic
// formula, so narrow local extrema are never missed by fixed sampling.
func hermiteAxisBounds64(p0, m0, p1, m1 float64) (lo, hi float64) {
	lo = math.Min(p0, p1)
	hi = math.Max(p0, p1)

	// Scale all four controls by a common positive factor before building the
	// polynomial. The curve scales linearly and derivative roots stay the
	// same, but coefficients stay bounded, keeping the solve finite and well
	// conditioned for extremely small or huge finite inputs.
	s := math.Max(math.Abs(p0), math.Max(math.Abs(m0), math.Max(math.Abs(p1), math.Abs(m1))))
	if s == 0 {
		return lo, hi
	}
	p0 /= s
	m0 /= s
	p1 /= s
	m1 /= s

	// p(t)/s = a*t^3 + b*t^2 + c*t + d
	a := 2*p0 + m0 + m1 - 2*p1
	b := -3*p0 - 2*m0 - m1 + 3*p1
	c := m0
	eval := func(t float64) float64 { return s * (((a*t+b)*t+c)*t + p0) }

	consider := func(t float64) {
		if t > 0 && t < 1 {
			v := eval(t)
			lo = math.Min(lo, v)
			hi = math.Max(hi, v)
		}
	}

	// p'(t)/s = A*t^2 + B*t + C
	A := 3 * a
	B := 2 * b
	C := c

	if A == 0 {
		if B != 0 {
			consider(-C / B)
		}
		return lo, hi
	}

	disc := B*B - 4*A*C
	if disc < 0 {
		// A slightly negative discriminant can be a cancellation artifact at
		// a (nearly) repeated root. Snap it to zero relative to the magnitude
		// of its terms. Considering the resulting vertex is always safe: it is
		// a real curve point that can only move a bound towards the true
		// extremum, never beyond the global min/max.
		termScale := math.Abs(B*B) + math.Abs(4*A*C)
		if disc < -1e-12*(1+termScale) {
			return lo, hi
		}
		disc = 0
	}

	sqrtDisc := math.Sqrt(disc)
	// Numerically stable quadratic roots (Numerical Recipes):
	// q = -1/2 * (B + copysign(sqrtDisc, B)); roots q/A and C/q.
	q := -0.5 * (B + math.Copysign(sqrtDisc, B))
	if q != 0 {
		consider(q / A)
		consider(C / q)
	} else {
		consider(-B / (2 * A))
	}
	return lo, hi
}
