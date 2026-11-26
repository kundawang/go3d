package hermit

import (
	"math"
	"math/rand"
	"testing"

	"github.com/ungerik/go3d/float64/vec3"
)

func assertBoxNear(t *testing.T, name string, got vec3.Box, wantMin, wantMax vec3.T, tol float64) {
	t.Helper()
	for axis := 0; axis < 3; axis++ {
		if d := math.Abs(got.Min[axis] - wantMin[axis]); d > tol {
			t.Errorf("%s Min[%d] = %.12g, want %.12g (diff %.3g > tol %.3g)",
				name, axis, got.Min[axis], wantMin[axis], d, tol)
		}
		if d := math.Abs(got.Max[axis] - wantMax[axis]); d > tol {
			t.Errorf("%s Max[%d] = %.12g, want %.12g (diff %.3g > tol %.3g)",
				name, axis, got.Max[axis], wantMax[axis], d, tol)
		}
		if got.Min[axis] > got.Max[axis] {
			t.Errorf("%s inverted box on axis %d: min %v > max %v", name, axis, got.Min[axis], got.Max[axis])
		}
	}
}

func assertPanics(t *testing.T, name string, f func()) {
	t.Helper()
	defer func() {
		if recover() == nil {
			t.Errorf("%s: expected panic, got none", name)
		}
	}()
	f()
}

// The motivating example: both endpoints at the origin and tangents
// (4,0,0) / (-4,0,0) give x(t) = 4*t*(1-t), range [0,1].
func TestBoundingBoxCanonicalArch64(t *testing.T) {
	herm := T{
		A: PointTangent{Point: vec3.T{0, 0, 0}, Tangent: vec3.T{4, 0, 0}},
		B: PointTangent{Point: vec3.T{0, 0, 0}, Tangent: vec3.T{-4, 0, 0}},
	}

	box := herm.BoundingBox()

	if box.Min[0] != 0 || box.Max[0] != 1 {
		t.Errorf("canonical X range = [%v, %v], want exactly [0,1]", box.Min[0], box.Max[0])
	}
	if box.Min[1] != 0 || box.Max[1] != 0 || box.Min[2] != 0 || box.Max[2] != 0 {
		t.Errorf("canonical Y/Z should be [0,0], got box %v", box)
	}

	// All entry points must agree.
	if m := herm.Bounds(); m != box {
		t.Errorf("Bounds() = %v, BoundingBox() = %v", m, box)
	}
	if m := BoundingBox(&herm.A.Point, &herm.A.Tangent, &herm.B.Point, &herm.B.Tangent); m != box {
		t.Errorf("package BoundingBox = %v, method = %v", m, box)
	}
	if m := Bounds(&herm.A.Point, &herm.A.Tangent, &herm.B.Point, &herm.B.Tangent); m != box {
		t.Errorf("package Bounds = %v, method = %v", m, box)
	}
	// Direct call on a composite literal must compile (value receiver).
	if m := (T{A: herm.A, B: herm.B}).BoundingBox(); m != box {
		t.Errorf("literal BoundingBox = %v, want %v", m, box)
	}
}

// Elevated Y/Z tracks: X is the canonical arch, Y is a downward quadratic
// with vertex 2 at t=0.25 and endpoint -2.5, Z is linear from -7 to -5.
func TestBoundingBoxElevatedYZ64(t *testing.T) {
	herm := T{
		A: PointTangent{Point: vec3.T{0, 1.5, -7}, Tangent: vec3.T{4, 4, 2}},
		B: PointTangent{Point: vec3.T{0, -2.5, -5}, Tangent: vec3.T{-4, -12, 2}},
	}
	box := herm.BoundingBox()
	wantMin := vec3.T{0, -2.5, -7}
	wantMax := vec3.T{1, 2, -5}
	assertBoxNear(t, "elevated", box, wantMin, wantMax, 1e-12)
}

// f(t) = 6t^3 - 9t^2 + 3t has both endpoints 0 and both tangents 3,
// but two interior extrema at (3 +/- sqrt(3))/6 with values +/-1/(2 sqrt(3)).
func TestBoundingBoxTwoInteriorExtrema64(t *testing.T) {
	herm := T{
		A: PointTangent{Point: vec3.T{0, 0, 0}, Tangent: vec3.T{3, 0, 0}},
		B: PointTangent{Point: vec3.T{0, 0, 0}, Tangent: vec3.T{3, 0, 0}},
	}
	box := herm.BoundingBox()
	want := 1.0 / (2.0 * math.Sqrt(3))
	assertBoxNear(t, "cubic", box, vec3.T{-want, 0, 0}, vec3.T{want, 0, 0}, 1e-12)
}

// A narrow peak at t=0.37 reaches 1 while a coarse fixed grid only sees 0.96.
func TestBoundingBoxNarrowPeakBeatsFixedSampling64(t *testing.T) {
	const r = 0.37
	const K = 100.0
	// f(t) = 1 - K*(t-r)^2 = -K t^2 + 2 K r t + (1 - K r^2)
	b := -K
	c := 2 * K * r
	d := 1 - K*r*r
	pA := d
	mA := c
	pB := -K + c + d
	mB := 2*b + c
	herm := T{
		A: PointTangent{Point: vec3.T{pA, 0, 0}, Tangent: vec3.T{mA, 0, 0}},
		B: PointTangent{Point: vec3.T{pB, 0, 0}, Tangent: vec3.T{mB, 0, 0}},
	}
	box := herm.BoundingBox()
	if math.Abs(box.Max[0]-1) > 1e-9 {
		t.Fatalf("narrow peak Max.X = %.12g, want ~1", box.Max[0])
	}

	coarseMax := -math.MaxFloat64
	for i := 0; i <= 20; i++ {
		p := herm.Point(float64(i) / 20)
		coarseMax = math.Max(coarseMax, p[0])
	}
	if coarseMax >= 0.97 {
		t.Fatalf("test setup: coarse grid unexpectedly caught peak: %v", coarseMax)
	}
	if box.Max[0] <= coarseMax+0.02 {
		t.Fatalf("analytic box %v not tighter than coarse sample max %v", box.Max[0], coarseMax)
	}
}

func TestBoundingBoxCollinear64(t *testing.T) {
	herm := T{
		A: PointTangent{Point: vec3.T{0, 0, 0}, Tangent: vec3.T{3, 4, 5}},
		B: PointTangent{Point: vec3.T{3, 4, 5}, Tangent: vec3.T{3, 4, 5}},
	}
	box := herm.BoundingBox()
	assertBoxNear(t, "line", box, herm.A.Point, herm.B.Point, 1e-12)
}

func TestBoundingBoxDegenerate64(t *testing.T) {
	zero := T{
		A: PointTangent{Point: vec3.T{}, Tangent: vec3.T{}},
		B: PointTangent{Point: vec3.T{}, Tangent: vec3.T{}},
	}
	box := zero.BoundingBox()
	if box.Min != (vec3.T{}) || box.Max != (vec3.T{}) {
		t.Errorf("zero curve box = %v, want zero box", box)
	}

	// Repeated point with zero tangents at a non-origin location.
	at := T{
		A: PointTangent{Point: vec3.T{-7.5, 2.5, 1}, Tangent: vec3.T{}},
		B: PointTangent{Point: vec3.T{-7.5, 2.5, 1}, Tangent: vec3.T{}},
	}
	b2 := at.BoundingBox()
	if b2.Min != at.A.Point || b2.Max != at.A.Point {
		t.Errorf("point curve box = %v, want %v", b2, at.A.Point)
	}
}

func TestBoundingBoxNegativeCoordinates64(t *testing.T) {
	herm := T{
		A: PointTangent{Point: vec3.T{-4, -4, -4}, Tangent: vec3.T{-8, 2, 0}},
		B: PointTangent{Point: vec3.T{-1, -9, -2}, Tangent: vec3.T{4, -2, 3}},
	}
	box := herm.BoundingBox()
	// Endpoints must be contained, box must never invert.
	for axis := 0; axis < 3; axis++ {
		if box.Min[axis] > herm.A.Point[axis]+1e-12 || box.Min[axis] > herm.B.Point[axis]+1e-12 {
			t.Errorf("Min[%d]=%v does not cover endpoints %v %v", axis, box.Min[axis], herm.A.Point[axis], herm.B.Point[axis])
		}
		if box.Max[axis] < herm.A.Point[axis]-1e-12 || box.Max[axis] < herm.B.Point[axis]-1e-12 {
			t.Errorf("Max[%d]=%v does not cover endpoints %v %v", axis, box.Max[axis], herm.A.Point[axis], herm.B.Point[axis])
		}
	}
}

func TestBoundingBoxTinyScale64(t *testing.T) {
	for _, scale := range []float64{1e-300, 1e-100, 1e-12, 1, 1e100, 1e300} {
		herm := T{
			A: PointTangent{Point: vec3.T{}, Tangent: vec3.T{4 * scale, 0, 0}},
			B: PointTangent{Point: vec3.T{}, Tangent: vec3.T{-4 * scale, 0, 0}},
		}
		box, err := herm.TryBoundingBox()
		if err != nil {
			t.Fatalf("scale %v: unexpected error %v", scale, err)
		}
		if math.Abs(box.Max[0]-scale) > 1e-9*scale {
			t.Errorf("scale %v: Max.X = %.12g, want %v", scale, box.Max[0], scale)
		}
		if box.Min[0] != 0 {
			t.Errorf("scale %v: Min.X = %.12g, want 0", scale, box.Min[0])
		}
	}
}

func TestBoundingBoxNonFinite64(t *testing.T) {
	bad := []vec3.T{
		{math.NaN(), 0, 0},
		{math.Inf(1), 0, 0},
		{0, math.Inf(-1), 0},
	}
	for i, v := range bad {
		herm := T{
			A: PointTangent{Point: v, Tangent: vec3.T{1, 1, 1}},
			B: PointTangent{Point: vec3.T{1, 1, 1}, Tangent: vec3.T{1, 1, 1}},
		}
		if _, err := herm.TryBoundingBox(); err == nil {
			t.Errorf("case %d: TryBoundingBox expected error for %v", i, v)
		}
		assertPanics(t, "BoundingBox non-finite", func() { herm.BoundingBox() })
	}

	if _, err := TryBoundingBox(nil, &vec3.T{}, &vec3.T{}, &vec3.T{}); err == nil {
		t.Error("nil pointA: expected error")
	}
	assertPanics(t, "package BoundingBox error", func() {
		BoundingBox(&vec3.T{math.NaN(), 0, 0}, &vec3.T{}, &vec3.T{}, &vec3.T{})
	})
}

// oracleAxisBox computes exact float64 polynomial extrema from independent
// coefficients and derivative roots.
func oracleAxisBox(p0, m0, p1, m1 float64) (lo, hi float64) {
	a := 2*p0 + m0 + m1 - 2*p1
	b := -3*p0 - 2*m0 - m1 + 3*p1
	c := m0
	eval := func(t float64) float64 { return ((a*t+b)*t+c)*t + p0 }
	lo = math.Min(p0, p1)
	hi = math.Max(p0, p1)
	consider := func(r float64) {
		if r > 0 && r < 1 {
			v := eval(r)
			lo = math.Min(lo, v)
			hi = math.Max(hi, v)
		}
	}
	A, B, C := 3*a, 2*b, c
	if A == 0 {
		if B != 0 {
			consider(-C / B)
		}
		return
	}
	disc := B*B - 4*A*C
	if disc >= 0 {
		s := math.Sqrt(disc)
		consider((-B - s) / (2 * A))
		consider((-B + s) / (2 * A))
	}
	return
}

func TestBoundingBoxRandomized64(t *testing.T) {
	rng := rand.New(rand.NewSource(20261127))
	scales := []float64{1e-30, 1e-8, 1, 1e8, 1e20}

	for trial := 0; trial < 400; trial++ {
		scale := scales[rng.Intn(len(scales))]
		rv := func() float64 { return scale * (2*rng.Float64() - 1) }
		herm := T{
			A: PointTangent{
				Point:   vec3.T{rv(), rv(), rv()},
				Tangent: vec3.T{rv(), rv(), rv()},
			},
			B: PointTangent{
				Point:   vec3.T{rv(), rv(), rv()},
				Tangent: vec3.T{rv(), rv(), rv()},
			},
		}
		// Every few trials force a degenerate pattern.
		switch trial % 7 {
		case 0:
			herm.B.Point = herm.A.Point
		case 1:
			herm.A.Tangent = vec3.T{}
			herm.B.Tangent = vec3.T{}
		case 2:
			herm.B.Point = herm.A.Point
			herm.A.Tangent = herm.B.Tangent
		}

		box, err := herm.TryBoundingBox()
		if err != nil {
			t.Fatalf("trial %d: unexpected error: %v", trial, err)
		}

		// Tightness against the analytic oracle.
		tol := 1e-8 * scale
		if scale == 0 {
			tol = 0
		}
		for axis := 0; axis < 3; axis++ {
			lo, hi := oracleAxisBox(
				herm.A.Point[axis], herm.A.Tangent[axis],
				herm.B.Point[axis], herm.B.Tangent[axis],
			)
			if d := math.Abs(box.Min[axis] - lo); d > tol {
				t.Fatalf("trial %d axis %d: Min %.12g want oracle %.12g (diff %g, tol %g)\nherm=%v",
					trial, axis, box.Min[axis], lo, d, tol, herm)
			}
			if d := math.Abs(box.Max[axis] - hi); d > tol {
				t.Fatalf("trial %d axis %d: Max %.12g want oracle %.12g (diff %g, tol %g)\nherm=%v",
					trial, axis, box.Max[axis], hi, d, tol, herm)
			}
		}

		// Containment of dense uniform + pseudo random float64 samples.
		cTol := 1e-9 * scale
		check := func(tt float64) {
			p := herm.Point(tt)
			for axis := 0; axis < 3; axis++ {
				if p[axis] < box.Min[axis]-cTol || p[axis] > box.Max[axis]+cTol {
					t.Fatalf("trial %d: Point(%v)[%d]=%.12g outside box [%g,%g] tol %g\nherm=%v box=%v",
						trial, tt, axis, p[axis], box.Min[axis], box.Max[axis], cTol, herm, box)
				}
			}
		}
		for i := 0; i <= 500; i++ {
			check(float64(i) / 500)
		}
		for i := 0; i < 500; i++ {
			check(rng.Float64())
		}
	}
}
