package hermit

import (
	"math"
	"testing"

	"github.com/ungerik/go3d/vec3"
)

const EPSILON = 0.0001

func TestParseAndString(t *testing.T) {
	original := T{
		A: PointTangent{
			Point:   vec3.T{1.0, 2.0, 3.0},
			Tangent: vec3.T{0.5, 0.5, 0.5},
		},
		B: PointTangent{
			Point:   vec3.T{4.0, 5.0, 6.0},
			Tangent: vec3.T{-0.5, 0.5, -0.5},
		},
	}

	// Convert to string and parse back
	str := original.String()
	parsed, err := Parse(str)
	if err != nil {
		t.Errorf("Parse() failed: %v", err)
	}

	// Check all components
	if math.Abs(float64(parsed.A.Point[0]-original.A.Point[0])) > EPSILON ||
		math.Abs(float64(parsed.A.Point[1]-original.A.Point[1])) > EPSILON ||
		math.Abs(float64(parsed.A.Point[2]-original.A.Point[2])) > EPSILON ||
		math.Abs(float64(parsed.A.Tangent[0]-original.A.Tangent[0])) > EPSILON ||
		math.Abs(float64(parsed.A.Tangent[1]-original.A.Tangent[1])) > EPSILON ||
		math.Abs(float64(parsed.A.Tangent[2]-original.A.Tangent[2])) > EPSILON ||
		math.Abs(float64(parsed.B.Point[0]-original.B.Point[0])) > EPSILON ||
		math.Abs(float64(parsed.B.Point[1]-original.B.Point[1])) > EPSILON ||
		math.Abs(float64(parsed.B.Point[2]-original.B.Point[2])) > EPSILON ||
		math.Abs(float64(parsed.B.Tangent[0]-original.B.Tangent[0])) > EPSILON ||
		math.Abs(float64(parsed.B.Tangent[1]-original.B.Tangent[1])) > EPSILON ||
		math.Abs(float64(parsed.B.Tangent[2]-original.B.Tangent[2])) > EPSILON {
		t.Errorf("Parse/String roundtrip failed: got %v, want %v", parsed, original)
	}
}

func TestPointEndpoints(t *testing.T) {
	// Create a simple hermit spline from (0,0,0) to (1,1,1) with tangents pointing in the same direction
	herm := T{
		A: PointTangent{
			Point:   vec3.T{0.0, 0.0, 0.0},
			Tangent: vec3.T{1.0, 1.0, 1.0},
		},
		B: PointTangent{
			Point:   vec3.T{1.0, 1.0, 1.0},
			Tangent: vec3.T{1.0, 1.0, 1.0},
		},
	}

	// At t=0, should return pointA
	p0 := herm.Point(0.0)
	if math.Abs(float64(p0[0]-herm.A.Point[0])) > EPSILON ||
		math.Abs(float64(p0[1]-herm.A.Point[1])) > EPSILON ||
		math.Abs(float64(p0[2]-herm.A.Point[2])) > EPSILON {
		t.Errorf("Point(0) should equal A.Point: got %v, want %v", p0, herm.A.Point)
	}

	// At t=1, should return pointB
	p1 := herm.Point(1.0)
	if math.Abs(float64(p1[0]-herm.B.Point[0])) > EPSILON ||
		math.Abs(float64(p1[1]-herm.B.Point[1])) > EPSILON ||
		math.Abs(float64(p1[2]-herm.B.Point[2])) > EPSILON {
		t.Errorf("Point(1) should equal B.Point: got %v, want %v", p1, herm.B.Point)
	}
}

func TestPointMidpoint(t *testing.T) {
	// Line along X-axis from (0,0,0) to (2,0,0) with X-aligned tangents
	herm := T{
		A: PointTangent{
			Point:   vec3.T{0.0, 0.0, 0.0},
			Tangent: vec3.T{1.0, 0.0, 0.0},
		},
		B: PointTangent{
			Point:   vec3.T{2.0, 0.0, 0.0},
			Tangent: vec3.T{1.0, 0.0, 0.0},
		},
	}

	// At t=0.5, should be approximately at (1, 0, 0) for a straight line
	p := herm.Point(0.5)

	// The y and z coordinates should be exactly 0
	if math.Abs(float64(p[1])) > EPSILON || math.Abs(float64(p[2])) > EPSILON {
		t.Errorf("Point(0.5) y,z coordinates should be 0: got %v", p)
	}

	// The x-coordinate should be approximately 1.0
	if math.Abs(float64(p[0]-1.0)) > EPSILON {
		t.Errorf("Point(0.5) x-coordinate should be ~1.0: got %v", p[0])
	}
}

func TestPointFunctionMatchesMethod(t *testing.T) {
	herm := T{
		A: PointTangent{
			Point:   vec3.T{0.0, 0.0, 0.0},
			Tangent: vec3.T{1.0, 0.5, 0.2},
		},
		B: PointTangent{
			Point:   vec3.T{2.0, 1.0, 0.5},
			Tangent: vec3.T{1.0, -0.5, 0.3},
		},
	}

	tests := []float32{0.0, 0.25, 0.5, 0.75, 1.0}
	for _, tVal := range tests {
		// Method call
		p1 := herm.Point(tVal)

		// Function call
		p2 := Point(&herm.A.Point, &herm.A.Tangent, &herm.B.Point, &herm.B.Tangent, tVal)

		if math.Abs(float64(p1[0]-p2[0])) > EPSILON ||
			math.Abs(float64(p1[1]-p2[1])) > EPSILON ||
			math.Abs(float64(p1[2]-p2[2])) > EPSILON {
			t.Errorf("Point method and function differ at t=%f: method=%v, func=%v", tVal, p1, p2)
		}
	}
}

func TestTangentAtZero(t *testing.T) {
	herm := T{
		A: PointTangent{
			Point:   vec3.T{0.0, 0.0, 0.0},
			Tangent: vec3.T{1.0, 0.5, 0.2},
		},
		B: PointTangent{
			Point:   vec3.T{2.0, 1.0, 0.5},
			Tangent: vec3.T{0.5, 1.0, 0.3},
		},
	}

	// At t=0, the Tangent function returns tangentA (verified by examining the formula at t=0)
	tan0 := herm.Tangent(0.0)
	expectedTan0 := herm.A.Tangent
	if math.Abs(float64(tan0[0]-expectedTan0[0])) > EPSILON ||
		math.Abs(float64(tan0[1]-expectedTan0[1])) > EPSILON ||
		math.Abs(float64(tan0[2]-expectedTan0[2])) > EPSILON {
		t.Errorf("Tangent(0) should equal A.Tangent: got %v, want %v", tan0, expectedTan0)
	}
}

func TestTangentFunctionMatchesMethod(t *testing.T) {
	herm := T{
		A: PointTangent{
			Point:   vec3.T{0.0, 0.0, 0.0},
			Tangent: vec3.T{1.0, 0.5, 0.2},
		},
		B: PointTangent{
			Point:   vec3.T{2.0, 1.0, 0.5},
			Tangent: vec3.T{0.5, 1.0, 0.3},
		},
	}

	tests := []float32{0.0, 0.25, 0.5, 0.75, 1.0}
	for _, tVal := range tests {
		// Method call
		tan1 := herm.Tangent(tVal)

		// Function call
		tan2 := Tangent(&herm.A.Point, &herm.A.Tangent, &herm.B.Point, &herm.B.Tangent, tVal)

		if math.Abs(float64(tan1[0]-tan2[0])) > EPSILON ||
			math.Abs(float64(tan1[1]-tan2[1])) > EPSILON ||
			math.Abs(float64(tan1[2]-tan2[2])) > EPSILON {
			t.Errorf("Tangent method and function differ at t=%f: method=%v, func=%v", tVal, tan1, tan2)
		}
	}
}

func TestLengthValues(t *testing.T) {
	herm := T{
		A: PointTangent{
			Point:   vec3.T{0.0, 0.0, 0.0},
			Tangent: vec3.T{1.0, 0.0, 0.0},
		},
		B: PointTangent{
			Point:   vec3.T{2.0, 0.0, 0.0},
			Tangent: vec3.T{1.0, 0.0, 0.0},
		},
	}

	// Test that Length returns a valid positive number
	length0 := herm.Length(0.0)
	if math.IsNaN(float64(length0)) || math.IsInf(float64(length0), 0) {
		t.Errorf("Length(0) should be a valid number: got %f", length0)
	}

	length1 := herm.Length(1.0)
	if math.IsNaN(float64(length1)) || math.IsInf(float64(length1), 0) {
		t.Errorf("Length(1) should be a valid number: got %f", length1)
	}
}

func TestLengthFunctionMatchesMethod(t *testing.T) {
	herm := T{
		A: PointTangent{
			Point:   vec3.T{0.0, 0.0, 0.0},
			Tangent: vec3.T{1.0, 0.5, 0.2},
		},
		B: PointTangent{
			Point:   vec3.T{2.0, 1.0, 0.5},
			Tangent: vec3.T{0.5, 1.0, 0.3},
		},
	}

	tests := []float32{0.0, 0.25, 0.5, 0.75, 1.0}
	for _, tVal := range tests {
		// Method call
		len1 := herm.Length(tVal)

		// Function call
		len2 := Length(&herm.A.Point, &herm.A.Tangent, &herm.B.Point, &herm.B.Tangent, tVal)

		if math.Abs(float64(len1-len2)) > EPSILON {
			t.Errorf("Length method and function differ at t=%f: method=%f, func=%f", tVal, len1, len2)
		}
	}
}

// straightLine returns the spline whose positions satisfy x(t) = 2t exactly:
// the control points are the axis endpoints and both tangents have magnitude 2,
// so the true arc length from 0 to t is 2t.
func straightLine() T {
	return T{
		A: PointTangent{
			Point:   vec3.T{0.0, 0.0, 0.0},
			Tangent: vec3.T{2.0, 0.0, 0.0},
		},
		B: PointTangent{
			Point:   vec3.T{2.0, 0.0, 0.0},
			Tangent: vec3.T{2.0, 0.0, 0.0},
		},
	}
}

// curved returns a genuinely non-straight spline used to check that the
// tangent is the true derivative and the length is the true arc length.
func curved() T {
	return T{
		A: PointTangent{
			Point:   vec3.T{0.0, 0.0, 0.0},
			Tangent: vec3.T{1.0, 2.0, -1.0},
		},
		B: PointTangent{
			Point:   vec3.T{3.0, 1.0, 2.0},
			Tangent: vec3.T{-1.0, 0.5, 2.0},
		},
	}
}

func TestLengthStraightLineExactValues(t *testing.T) {
	herm := straightLine()

	if got := herm.Length(0.5); math.Abs(float64(got-1.0)) > 1e-6 {
		t.Errorf("Length(0.5) = %v, want 1", got)
	}
	if got := herm.Length(1.0); math.Abs(float64(got-2.0)) > 1e-6 {
		t.Errorf("Length(1) = %v, want 2", got)
	}
	if got := herm.Length(0.0); got != 0 {
		t.Errorf("Length(0) = %v, want 0", got)
	}

	want := vec3.T{2.0, 0.0, 0.0}
	tan := herm.Tangent(0.5)
	for i := 0; i < 3; i++ {
		if math.Abs(float64(tan[i]-want[i])) > 1e-6 {
			t.Errorf("Tangent(0.5) = %v, want %v", tan, want)
			break
		}
	}
}

func TestTangentIsCurveDerivative(t *testing.T) {
	herm := curved()

	// Tangent at the endpoints must reproduce the given tangents exactly.
	if got := herm.Tangent(0.0); math.Abs(float64(got[0]-herm.A.Tangent[0])) > EPSILON ||
		math.Abs(float64(got[1]-herm.A.Tangent[1])) > EPSILON ||
		math.Abs(float64(got[2]-herm.A.Tangent[2])) > EPSILON {
		t.Errorf("Tangent(0) = %v, want A.Tangent %v", got, herm.A.Tangent)
	}
	if got := herm.Tangent(1.0); math.Abs(float64(got[0]-herm.B.Tangent[0])) > EPSILON ||
		math.Abs(float64(got[1]-herm.B.Tangent[1])) > EPSILON ||
		math.Abs(float64(got[2]-herm.B.Tangent[2])) > EPSILON {
		t.Errorf("Tangent(1) = %v, want B.Tangent %v", got, herm.B.Tangent)
	}

	// Compare against a central difference of Point for interior samples.
	const h = 1e-3
	for _, tv := range []float32{0.2, 0.45, 0.5, 0.75} {
		lo := herm.Point(tv - h)
		hi := herm.Point(tv + h)
		var numeric vec3.T
		for i := 0; i < 3; i++ {
			numeric[i] = (hi[i] - lo[i]) / (2 * h)
		}
		got := herm.Tangent(tv)
		for i := 0; i < 3; i++ {
			if math.Abs(float64(got[i]-numeric[i])) > 5e-3 {
				t.Errorf("Tangent(%v)[%d] = %v, numerical derivative %v", tv, i, got[i], numeric[i])
			}
		}
	}
}

func TestLengthNonStraight(t *testing.T) {
	herm := curved()

	// Independent reference: length of a densely sampled polyline through
	// Point, accumulated in float64.
	const samples = 200000
	ref := 0.0
	prev := herm.Point(0)
	for i := 1; i <= samples; i++ {
		cur := herm.Point(float32(float64(i) / samples))
		var seg vec3.T
		for d := 0; d < 3; d++ {
			seg[d] = cur[d] - prev[d]
		}
		ref += float64(seg.Length())
		prev = cur
	}

	if got := float64(herm.Length(1.0)); math.Abs(got-ref) > 1e-4*ref {
		t.Errorf("Length(1) = %v, reference arc length %v", got, ref)
	}

	// Anchoring at 0 and monotonicity across a fine sweep.
	last := 0.0
	for i := 0; i <= 1000; i++ {
		l := float64(herm.Length(float32(i) / 1000))
		if l < last {
			t.Fatalf("Length is not monotone at t=%v: %v < %v", float32(i)/1000, l, last)
		}
		last = l
	}
	if math.Abs(last-float64(herm.Length(1.0))) > 1e-6 {
		t.Errorf("fine sweep did not reach Length(1): %v vs %v", last, herm.Length(1.0))
	}
}

func TestLengthDegenerateCurve(t *testing.T) {
	var herm T // all control values zero, derivative identically zero

	if got := herm.Length(0.5); got != 0 {
		t.Errorf("degenerate Length(0.5) = %v, want 0", got)
	}
	if got := herm.Length(1.0); got != 0 {
		t.Errorf("degenerate Length(1) = %v, want 0", got)
	}
	tan := herm.Tangent(0.5)
	for i := 0; i < 3; i++ {
		if tan[i] != 0 {
			t.Errorf("degenerate Tangent(0.5) = %v, want zero", tan)
			break
		}
	}
}

func TestLengthDomainAndNaN(t *testing.T) {
	herm := curved()
	total := herm.Length(1.0)

	if got := herm.Length(-0.5); got != 0 {
		t.Errorf("Length(-0.5) = %v, want 0", got)
	}
	if got := herm.Length(1.5); math.Abs(float64(got-total)) > 1e-6 {
		t.Errorf("Length(1.5) = %v, want clamped total %v", got, total)
	}
	if got := herm.Length(float32(math.NaN())); !math.IsNaN(float64(got)) {
		t.Errorf("Length(NaN) = %v, want NaN", got)
	}
	if got := herm.Tangent(float32(math.NaN())); !math.IsNaN(float64(got[0])) {
		t.Errorf("Tangent(NaN) = %v, want NaN", got)
	}

	// Extrapolation beyond the domain stays consistent with Point.
	const h = 1e-3
	lo := herm.Point(1.2 - h)
	hi := herm.Point(1.2 + h)
	numeric := (hi[0] - lo[0]) / (2 * h)
	if got := herm.Tangent(1.2)[0]; math.Abs(float64(got-numeric)) > 5e-3 {
		t.Errorf("Tangent(1.2)[0] = %v, numerical derivative %v", got, numeric)
	}
}

func TestNonFiniteControls(t *testing.T) {
	nan := float32(math.NaN())
	inf := float32(math.Inf(1))
	finite := straightLine()

	cases := []struct {
		name string
		make func() T
	}{
		{"NaN pointA", func() T { h := finite; h.A.Point[1] = nan; return h }},
		{"NaN tangentA", func() T { h := finite; h.A.Tangent[1] = nan; return h }},
		{"Inf tangentB", func() T { h := finite; h.B.Tangent[2] = inf; return h }},
		{"Inf pointB", func() T { h := finite; h.B.Point[0] = inf; return h }},
	}
	for _, c := range cases {
		herm := c.make()
		if got := herm.Length(0.5); !math.IsNaN(float64(got)) {
			t.Errorf("%s: Length(0.5) = %v, want NaN", c.name, got)
		}
		if got := herm.Tangent(0.5); !math.IsNaN(float64(got[0])) {
			t.Errorf("%s: Tangent(0.5) = %v, want NaN", c.name, got)
		}
	}
}
