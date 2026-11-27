package hermit2

import (
	"math"
	"testing"

	"github.com/ungerik/go3d/vec2"
)

const EPSILON = 0.0001

func TestParseAndString(t *testing.T) {
	original := T{
		A: PointTangent{
			Point:   vec2.T{1.0, 2.0},
			Tangent: vec2.T{0.5, 0.5},
		},
		B: PointTangent{
			Point:   vec2.T{3.0, 4.0},
			Tangent: vec2.T{-0.5, 0.5},
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
		math.Abs(float64(parsed.A.Tangent[0]-original.A.Tangent[0])) > EPSILON ||
		math.Abs(float64(parsed.A.Tangent[1]-original.A.Tangent[1])) > EPSILON ||
		math.Abs(float64(parsed.B.Point[0]-original.B.Point[0])) > EPSILON ||
		math.Abs(float64(parsed.B.Point[1]-original.B.Point[1])) > EPSILON ||
		math.Abs(float64(parsed.B.Tangent[0]-original.B.Tangent[0])) > EPSILON ||
		math.Abs(float64(parsed.B.Tangent[1]-original.B.Tangent[1])) > EPSILON {
		t.Errorf("Parse/String roundtrip failed: got %v, want %v", parsed, original)
	}
}

func TestPointEndpoints(t *testing.T) {
	// Create a simple hermit spline from (0,0) to (1,1) with tangents pointing in the same direction
	herm := T{
		A: PointTangent{
			Point:   vec2.T{0.0, 0.0},
			Tangent: vec2.T{1.0, 1.0},
		},
		B: PointTangent{
			Point:   vec2.T{1.0, 1.0},
			Tangent: vec2.T{1.0, 1.0},
		},
	}

	// At t=0, should return pointA
	p0 := herm.Point(0.0)
	if math.Abs(float64(p0[0]-herm.A.Point[0])) > EPSILON ||
		math.Abs(float64(p0[1]-herm.A.Point[1])) > EPSILON {
		t.Errorf("Point(0) should equal A.Point: got %v, want %v", p0, herm.A.Point)
	}

	// At t=1, should return pointB
	p1 := herm.Point(1.0)
	if math.Abs(float64(p1[0]-herm.B.Point[0])) > EPSILON ||
		math.Abs(float64(p1[1]-herm.B.Point[1])) > EPSILON {
		t.Errorf("Point(1) should equal B.Point: got %v, want %v", p1, herm.B.Point)
	}
}

func TestPointMidpoint(t *testing.T) {
	// Horizontal line from (0,0) to (2,0) with horizontal tangents
	herm := T{
		A: PointTangent{
			Point:   vec2.T{0.0, 0.0},
			Tangent: vec2.T{1.0, 0.0},
		},
		B: PointTangent{
			Point:   vec2.T{2.0, 0.0},
			Tangent: vec2.T{1.0, 0.0},
		},
	}

	// At t=0.5, should be approximately at (1, 0) for a straight line
	p := herm.Point(0.5)

	// The y-coordinate should be exactly 0 for horizontal tangents
	if math.Abs(float64(p[1])) > EPSILON {
		t.Errorf("Point(0.5) y-coordinate should be 0: got %v", p[1])
	}

	// The x-coordinate should be approximately 1.0
	if math.Abs(float64(p[0]-1.0)) > EPSILON {
		t.Errorf("Point(0.5) x-coordinate should be ~1.0: got %v", p[0])
	}
}

func TestPointFunctionMatchesMethod(t *testing.T) {
	herm := T{
		A: PointTangent{
			Point:   vec2.T{0.0, 0.0},
			Tangent: vec2.T{1.0, 0.5},
		},
		B: PointTangent{
			Point:   vec2.T{2.0, 1.0},
			Tangent: vec2.T{1.0, -0.5},
		},
	}

	tests := []float32{0.0, 0.25, 0.5, 0.75, 1.0}
	for _, tVal := range tests {
		// Method call
		p1 := herm.Point(tVal)

		// Function call
		p2 := Point(&herm.A.Point, &herm.A.Tangent, &herm.B.Point, &herm.B.Tangent, tVal)

		if math.Abs(float64(p1[0]-p2[0])) > EPSILON ||
			math.Abs(float64(p1[1]-p2[1])) > EPSILON {
			t.Errorf("Point method and function differ at t=%f: method=%v, func=%v", tVal, p1, p2)
		}
	}
}

func TestTangentAtZero(t *testing.T) {
	herm := T{
		A: PointTangent{
			Point:   vec2.T{0.0, 0.0},
			Tangent: vec2.T{1.0, 0.5},
		},
		B: PointTangent{
			Point:   vec2.T{2.0, 1.0},
			Tangent: vec2.T{0.5, 1.0},
		},
	}

	// At t=0, the Tangent function returns tangentA (verified by examining the formula at t=0)
	tan0 := herm.Tangent(0.0)
	expectedTan0 := herm.A.Tangent
	if math.Abs(float64(tan0[0]-expectedTan0[0])) > EPSILON ||
		math.Abs(float64(tan0[1]-expectedTan0[1])) > EPSILON {
		t.Errorf("Tangent(0) should equal A.Tangent: got %v, want %v", tan0, expectedTan0)
	}
}

func TestTangentFunctionMatchesMethod(t *testing.T) {
	herm := T{
		A: PointTangent{
			Point:   vec2.T{0.0, 0.0},
			Tangent: vec2.T{1.0, 0.5},
		},
		B: PointTangent{
			Point:   vec2.T{2.0, 1.0},
			Tangent: vec2.T{0.5, 1.0},
		},
	}

	tests := []float32{0.0, 0.25, 0.5, 0.75, 1.0}
	for _, tVal := range tests {
		// Method call
		tan1 := herm.Tangent(tVal)

		// Function call
		tan2 := Tangent(&herm.A.Point, &herm.A.Tangent, &herm.B.Point, &herm.B.Tangent, tVal)

		if math.Abs(float64(tan1[0]-tan2[0])) > EPSILON ||
			math.Abs(float64(tan1[1]-tan2[1])) > EPSILON {
			t.Errorf("Tangent method and function differ at t=%f: method=%v, func=%v", tVal, tan1, tan2)
		}
	}
}

func TestLengthValues(t *testing.T) {
	herm := T{
		A: PointTangent{
			Point:   vec2.T{0.0, 0.0},
			Tangent: vec2.T{1.0, 0.0},
		},
		B: PointTangent{
			Point:   vec2.T{2.0, 0.0},
			Tangent: vec2.T{1.0, 0.0},
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
			Point:   vec2.T{0.0, 0.0},
			Tangent: vec2.T{1.0, 0.5},
		},
		B: PointTangent{
			Point:   vec2.T{2.0, 1.0},
			Tangent: vec2.T{0.5, 1.0},
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

func straightLine2() T {
	return T{
		A: PointTangent{
			Point:   vec2.T{0.0, 0.0},
			Tangent: vec2.T{2.0, 0.0},
		},
		B: PointTangent{
			Point:   vec2.T{2.0, 0.0},
			Tangent: vec2.T{2.0, 0.0},
		},
	}
}

func curved2() T {
	return T{
		A: PointTangent{
			Point:   vec2.T{0.0, 0.0},
			Tangent: vec2.T{1.0, 2.0},
		},
		B: PointTangent{
			Point:   vec2.T{3.0, 1.0},
			Tangent: vec2.T{-1.0, 0.5},
		},
	}
}

func TestLengthStraightLineExactValues(t *testing.T) {
	herm := straightLine2()

	if got := herm.Length(0.5); math.Abs(float64(got-1.0)) > 1e-6 {
		t.Errorf("Length(0.5) = %v, want 1", got)
	}
	if got := herm.Length(1.0); math.Abs(float64(got-2.0)) > 1e-6 {
		t.Errorf("Length(1) = %v, want 2", got)
	}
	if got := herm.Length(0.0); got != 0 {
		t.Errorf("Length(0) = %v, want 0", got)
	}

	tan := herm.Tangent(0.5)
	if math.Abs(float64(tan[0]-2.0)) > 1e-6 || math.Abs(float64(tan[1])) > 1e-6 {
		t.Errorf("Tangent(0.5) = %v, want (2,0)", tan)
	}
}

func TestTangentIsCurveDerivative(t *testing.T) {
	herm := curved2()

	if got := herm.Tangent(0.0); math.Abs(float64(got[0]-herm.A.Tangent[0])) > EPSILON ||
		math.Abs(float64(got[1]-herm.A.Tangent[1])) > EPSILON {
		t.Errorf("Tangent(0) = %v, want A.Tangent %v", got, herm.A.Tangent)
	}
	if got := herm.Tangent(1.0); math.Abs(float64(got[0]-herm.B.Tangent[0])) > EPSILON ||
		math.Abs(float64(got[1]-herm.B.Tangent[1])) > EPSILON {
		t.Errorf("Tangent(1) = %v, want B.Tangent %v", got, herm.B.Tangent)
	}

	const h = 1e-3
	for _, tv := range []float32{0.2, 0.45, 0.5, 0.75} {
		lo := herm.Point(tv - h)
		hi := herm.Point(tv + h)
		var numeric vec2.T
		for i := 0; i < 2; i++ {
			numeric[i] = (hi[i] - lo[i]) / (2 * h)
		}
		got := herm.Tangent(tv)
		for i := 0; i < 2; i++ {
			if math.Abs(float64(got[i]-numeric[i])) > 5e-3 {
				t.Errorf("Tangent(%v)[%d] = %v, numerical derivative %v", tv, i, got[i], numeric[i])
			}
		}
	}
}

func TestLengthNonStraight(t *testing.T) {
	herm := curved2()

	const samples = 200000
	ref := 0.0
	prev := herm.Point(0)
	for i := 1; i <= samples; i++ {
		cur := herm.Point(float32(float64(i) / samples))
		var seg vec2.T
		for d := 0; d < 2; d++ {
			seg[d] = cur[d] - prev[d]
		}
		ref += float64(seg.Length())
		prev = cur
	}
	if got := float64(herm.Length(1.0)); math.Abs(got-ref) > 1e-4*ref {
		t.Errorf("Length(1) = %v, reference arc length %v", got, ref)
	}

	last := 0.0
	for i := 0; i <= 1000; i++ {
		l := float64(herm.Length(float32(i) / 1000))
		if l < last {
			t.Fatalf("Length is not monotone at t=%v: %v < %v", float32(i)/1000, l, last)
		}
		last = l
	}
}

func TestLengthDegenerateCurve(t *testing.T) {
	var herm T
	if got := herm.Length(0.5); got != 0 {
		t.Errorf("degenerate Length(0.5) = %v, want 0", got)
	}
	tan := herm.Tangent(0.5)
	if tan[0] != 0 || tan[1] != 0 {
		t.Errorf("degenerate Tangent(0.5) = %v, want zero", tan)
	}
}

func TestLengthDomainAndNaN(t *testing.T) {
	herm := curved2()
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
}

func TestNonFiniteControls(t *testing.T) {
	nan := float32(math.NaN())
	inf := float32(math.Inf(1))
	finite := straightLine2()

	cases := []struct {
		name string
		make func() T
	}{
		{"NaN pointA", func() T { h := finite; h.A.Point[1] = nan; return h }},
		{"NaN tangentA", func() T { h := finite; h.A.Tangent[1] = nan; return h }},
		{"Inf tangentB", func() T { h := finite; h.B.Tangent[1] = inf; return h }},
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
