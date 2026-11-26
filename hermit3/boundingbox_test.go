package hermit

import (
	"math"
	"math/rand"
	"testing"

	f64hermit "github.com/ungerik/go3d/float64/hermit3"
	f64vec3 "github.com/ungerik/go3d/float64/vec3"
	"github.com/ungerik/go3d/vec3"
)

func assertBoxNear32(t *testing.T, name string, got vec3.Box, wantMin, wantMax vec3.T, tol float32) {
	t.Helper()
	for axis := 0; axis < 3; axis++ {
		if d := float32(math.Abs(float64(got.Min[axis] - wantMin[axis]))); d > tol {
			t.Errorf("%s Min[%d] = %.10g, want %.10g (diff %.3g > tol %.3g)",
				name, axis, got.Min[axis], wantMin[axis], d, tol)
		}
		if d := float32(math.Abs(float64(got.Max[axis] - wantMax[axis]))); d > tol {
			t.Errorf("%s Max[%d] = %.10g, want %.10g (diff %.3g > tol %.3g)",
				name, axis, got.Max[axis], wantMax[axis], d, tol)
		}
		if got.Min[axis] > got.Max[axis] {
			t.Errorf("%s inverted box on axis %d: min %v > max %v", name, axis, got.Min[axis], got.Max[axis])
		}
	}
}

func assertPanics32(t *testing.T, name string, f func()) {
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
func TestBoundingBoxCanonicalArch32(t *testing.T) {
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

	if m := herm.Bounds(); m != box {
		t.Errorf("Bounds() = %v, BoundingBox() = %v", m, box)
	}
	if m := BoundingBox(&herm.A.Point, &herm.A.Tangent, &herm.B.Point, &herm.B.Tangent); m != box {
		t.Errorf("package BoundingBox = %v, method = %v", m, box)
	}
	if m := Bounds(&herm.A.Point, &herm.A.Tangent, &herm.B.Point, &herm.B.Tangent); m != box {
		t.Errorf("package Bounds = %v, method = %v", m, box)
	}
	if m := (T{A: herm.A, B: herm.B}).BoundingBox(); m != box {
		t.Errorf("literal BoundingBox = %v, want %v", m, box)
	}
}

// Elevated Y/Z tracks: X is the canonical arch, Y is a downward quadratic
// with vertex 2 at t=0.25 and endpoint -2.5, Z is linear from -7 to -5.
func TestBoundingBoxElevatedYZ32(t *testing.T) {
	herm := T{
		A: PointTangent{Point: vec3.T{0, 1.5, -7}, Tangent: vec3.T{4, 4, 2}},
		B: PointTangent{Point: vec3.T{0, -2.5, -5}, Tangent: vec3.T{-4, -12, 2}},
	}
	box := herm.BoundingBox()
	wantMin := vec3.T{0, -2.5, -7}
	wantMax := vec3.T{1, 2, -5}
	assertBoxNear32(t, "elevated", box, wantMin, wantMax, 1e-6)
}

// f(t) = 6t^3 - 9t^2 + 3t has both endpoints 0 and both tangents 3,
// but two interior extrema at (3 +/- sqrt(3))/6 with values +/-1/(2 sqrt(3)).
func TestBoundingBoxTwoInteriorExtrema32(t *testing.T) {
	herm := T{
		A: PointTangent{Point: vec3.T{0, 0, 0}, Tangent: vec3.T{3, 0, 0}},
		B: PointTangent{Point: vec3.T{0, 0, 0}, Tangent: vec3.T{3, 0, 0}},
	}
	box := herm.BoundingBox()
	want := float32(1.0 / (2.0 * math.Sqrt(3)))
	assertBoxNear32(t, "cubic", box, vec3.T{-want, 0, 0}, vec3.T{want, 0, 0}, 1e-6)
}

// A narrow peak at t=0.37 reaches 1 while a coarse fixed grid only sees 0.96.
func TestBoundingBoxNarrowPeakBeatsFixedSampling32(t *testing.T) {
	const r = 0.37
	const K = 100.0
	b := float32(-K)
	c := float32(2 * K * r)
	d := float32(1 - K*r*r)
	pA := d
	mA := c
	pB := b + c + d
	mB := 2*b + c
	herm := T{
		A: PointTangent{Point: vec3.T{pA, 0, 0}, Tangent: vec3.T{mA, 0, 0}},
		B: PointTangent{Point: vec3.T{pB, 0, 0}, Tangent: vec3.T{mB, 0, 0}},
	}
	box := herm.BoundingBox()
	if math.Abs(float64(box.Max[0])-1) > 1e-5 {
		t.Fatalf("narrow peak Max.X = %.10g, want ~1", box.Max[0])
	}

	var coarseMax float32 = -math.MaxFloat32
	for i := 0; i <= 20; i++ {
		p := herm.Point(float32(i) / 20)
		coarseMax = max32(coarseMax, p[0])
	}
	if coarseMax >= 0.97 {
		t.Fatalf("test setup: coarse grid unexpectedly caught peak: %v", coarseMax)
	}
	if box.Max[0] <= coarseMax+0.02 {
		t.Fatalf("analytic box %v not tighter than coarse sample max %v", box.Max[0], coarseMax)
	}
}

func max32(a, b float32) float32 {
	if a > b {
		return a
	}
	return b
}

func TestBoundingBoxCollinear32(t *testing.T) {
	herm := T{
		A: PointTangent{Point: vec3.T{0, 0, 0}, Tangent: vec3.T{3, 4, 5}},
		B: PointTangent{Point: vec3.T{3, 4, 5}, Tangent: vec3.T{3, 4, 5}},
	}
	box := herm.BoundingBox()
	assertBoxNear32(t, "line", box, herm.A.Point, herm.B.Point, 1e-6)
}

func TestBoundingBoxDegenerate32(t *testing.T) {
	zero := T{
		A: PointTangent{Point: vec3.T{}, Tangent: vec3.T{}},
		B: PointTangent{Point: vec3.T{}, Tangent: vec3.T{}},
	}
	box := zero.BoundingBox()
	if box.Min != (vec3.T{}) || box.Max != (vec3.T{}) {
		t.Errorf("zero curve box = %v, want zero box", box)
	}

	at := T{
		A: PointTangent{Point: vec3.T{-7.5, 2.5, 1}, Tangent: vec3.T{}},
		B: PointTangent{Point: vec3.T{-7.5, 2.5, 1}, Tangent: vec3.T{}},
	}
	b2 := at.BoundingBox()
	if b2.Min != at.A.Point || b2.Max != at.A.Point {
		t.Errorf("point curve box = %v, want %v", b2, at.A.Point)
	}
}

func TestBoundingBoxNegativeCoordinates32(t *testing.T) {
	herm := T{
		A: PointTangent{Point: vec3.T{-4, -4, -4}, Tangent: vec3.T{-8, 2, 0}},
		B: PointTangent{Point: vec3.T{-1, -9, -2}, Tangent: vec3.T{4, -2, 3}},
	}
	box := herm.BoundingBox()
	for axis := 0; axis < 3; axis++ {
		if box.Min[axis] > herm.A.Point[axis]+1e-6 || box.Min[axis] > herm.B.Point[axis]+1e-6 {
			t.Errorf("Min[%d]=%v does not cover endpoints %v %v", axis, box.Min[axis], herm.A.Point[axis], herm.B.Point[axis])
		}
		if box.Max[axis] < herm.A.Point[axis]-1e-6 || box.Max[axis] < herm.B.Point[axis]-1e-6 {
			t.Errorf("Max[%d]=%v does not cover endpoints %v %v", axis, box.Max[axis], herm.A.Point[axis], herm.B.Point[axis])
		}
	}
}

func TestBoundingBoxTinyScale32(t *testing.T) {
	for _, scale := range []float32{1e-30, 1e-20, 1e-12, 1, 1e20, 1e30} {
		herm := T{
			A: PointTangent{Point: vec3.T{}, Tangent: vec3.T{4 * scale, 0, 0}},
			B: PointTangent{Point: vec3.T{}, Tangent: vec3.T{-4 * scale, 0, 0}},
		}
		box, err := herm.TryBoundingBox()
		if err != nil {
			t.Fatalf("scale %v: unexpected error %v", scale, err)
		}
		if d := float32(math.Abs(float64(box.Max[0] - scale))); d > 1e-5*scale {
			t.Errorf("scale %v: Max.X = %.10g, want %v (rel tol)", scale, box.Max[0], scale)
		}
		if box.Min[0] != 0 {
			t.Errorf("scale %v: Min.X = %.10g, want 0", scale, box.Min[0])
		}
	}
}

func TestBoundingBoxNonFinite32(t *testing.T) {
	bad := []vec3.T{
		{float32(math.NaN()), 0, 0},
		{float32(math.Inf(1)), 0, 0},
		{0, float32(math.Inf(-1)), 0},
	}
	for i, v := range bad {
		herm := T{
			A: PointTangent{Point: v, Tangent: vec3.T{1, 1, 1}},
			B: PointTangent{Point: vec3.T{1, 1, 1}, Tangent: vec3.T{1, 1, 1}},
		}
		if _, err := herm.TryBoundingBox(); err == nil {
			t.Errorf("case %d: TryBoundingBox expected error for %v", i, v)
		}
		assertPanics32(t, "BoundingBox non-finite", func() { herm.BoundingBox() })
	}

	if _, err := TryBoundingBox(nil, &vec3.T{}, &vec3.T{}, &vec3.T{}); err == nil {
		t.Error("nil pointA: expected error")
	}
	assertPanics32(t, "package BoundingBox error", func() {
		BoundingBox(&vec3.T{float32(math.NaN()), 0, 0}, &vec3.T{}, &vec3.T{}, &vec3.T{})
	})
}

// Independent float64 analytic oracle of the exact float32-stored controls.
func oracleAxisBox32(p0, m0, p1, m1 float32) (lo, hi float64) {
	P0, M0, P1, M1 := float64(p0), float64(m0), float64(p1), float64(m1)
	a := 2*P0 + M0 + M1 - 2*P1
	b := -3*P0 - 2*M0 - M1 + 3*P1
	c := M0
	eval := func(t float64) float64 { return ((a*t+b)*t+c)*t + P0 }
	lo = math.Min(P0, P1)
	hi = math.Max(P0, P1)
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

func TestBoundingBoxRandomized32(t *testing.T) {
	rng := rand.New(rand.NewSource(20261127))
	scales := []float32{1e-20, 1e-8, 1, 1e8, 1e20}

	for trial := 0; trial < 600; trial++ {
		scale := scales[rng.Intn(len(scales))]
		rv := func() float32 { return scale * float32(2*rng.Float64()-1) }
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

		tol := 1e-5 * float64(scale)
		for axis := 0; axis < 3; axis++ {
			lo, hi := oracleAxisBox32(
				herm.A.Point[axis], herm.A.Tangent[axis],
				herm.B.Point[axis], herm.B.Tangent[axis],
			)
			if d := math.Abs(float64(box.Min[axis]) - lo); d > tol {
				t.Fatalf("trial %d axis %d: Min %.10g want oracle %.10g (diff %g, tol %g)\nherm=%v",
					trial, axis, box.Min[axis], lo, d, tol, herm)
			}
			if d := math.Abs(float64(box.Max[axis]) - hi); d > tol {
				t.Fatalf("trial %d axis %d: Max %.10g want oracle %.10g (diff %g, tol %g)\nherm=%v",
					trial, axis, box.Max[axis], hi, d, tol, herm)
			}
		}

		cTol := 1e-4 * scale
		check := func(tt float32) {
			p := herm.Point(tt)
			for axis := 0; axis < 3; axis++ {
				if p[axis] < box.Min[axis]-cTol || p[axis] > box.Max[axis]+cTol {
					t.Fatalf("trial %d: Point(%v)[%d]=%.10g outside box [%g,%g] tol %g\nherm=%v box=%v",
						trial, tt, axis, p[axis], box.Min[axis], box.Max[axis], cTol, herm, box)
				}
			}
		}
		for i := 0; i <= 500; i++ {
			check(float32(i) / 500)
		}
		for i := 0; i < 500; i++ {
			check(float32(rng.Float64()))
		}
	}
}

// The float32 and float64 packages must give corresponding boxes when the
// float64 side uses the exact float32-stored control values.
func TestBoundingBoxFloat64Correspondence32(t *testing.T) {
	cases := []T{
		{
			A: PointTangent{Point: vec3.T{}, Tangent: vec3.T{4, 0, 0}},
			B: PointTangent{Point: vec3.T{}, Tangent: vec3.T{-4, 0, 0}},
		},
		{
			A: PointTangent{Point: vec3.T{0, 1.5, -7}, Tangent: vec3.T{4, 4, 2}},
			B: PointTangent{Point: vec3.T{0, -2.5, -5}, Tangent: vec3.T{-4, -12, 2}},
		},
		{
			A: PointTangent{Point: vec3.T{0, 0, 0}, Tangent: vec3.T{3, 0, 0}},
			B: PointTangent{Point: vec3.T{0, 0, 0}, Tangent: vec3.T{3, 0, 0}},
		},
	}
	for ci, h32 := range cases {
		h64 := f64hermit.T{
			A: f64hermit.PointTangent{
				Point:   f64vec3.T{float64(h32.A.Point[0]), float64(h32.A.Point[1]), float64(h32.A.Point[2])},
				Tangent: f64vec3.T{float64(h32.A.Tangent[0]), float64(h32.A.Tangent[1]), float64(h32.A.Tangent[2])},
			},
			B: f64hermit.PointTangent{
				Point:   f64vec3.T{float64(h32.B.Point[0]), float64(h32.B.Point[1]), float64(h32.B.Point[2])},
				Tangent: f64vec3.T{float64(h32.B.Tangent[0]), float64(h32.B.Tangent[1]), float64(h32.B.Tangent[2])},
			},
		}
		b32 := h32.BoundingBox()
		b64 := h64.BoundingBox()
		for axis := 0; axis < 3; axis++ {
			tol := 1e-6 * math.Max(1, math.Abs(b64.Min[axis])+math.Abs(b64.Max[axis]))
			if d := math.Abs(float64(b32.Min[axis]) - b64.Min[axis]); d > tol {
				t.Errorf("case %d Min[%d]: float32 %v vs float64 %v diff %g", ci, axis, b32.Min[axis], b64.Min[axis], d)
			}
			if d := math.Abs(float64(b32.Max[axis]) - b64.Max[axis]); d > tol {
				t.Errorf("case %d Max[%d]: float32 %v vs float64 %v diff %g", ci, axis, b32.Max[axis], b64.Max[axis], d)
			}
		}
	}

	// Randomized correspondence over wide scales. The float64 side is built
	// from the exact float32-stored values.
	rng := rand.New(rand.NewSource(777))
	for trial := 0; trial < 300; trial++ {
		scale := []float32{1e-20, 1e-8, 1, 1e8, 1e20}[rng.Intn(5)]
		rv := func() float32 { return scale * float32(2*rng.Float64()-1) }
		h32 := T{
			A: PointTangent{Point: vec3.T{rv(), rv(), rv()}, Tangent: vec3.T{rv(), rv(), rv()}},
			B: PointTangent{Point: vec3.T{rv(), rv(), rv()}, Tangent: vec3.T{rv(), rv(), rv()}},
		}
		cv := func(v vec3.T) f64vec3.T {
			return f64vec3.T{float64(v[0]), float64(v[1]), float64(v[2])}
		}
		h64 := f64hermit.T{
			A: f64hermit.PointTangent{Point: cv(h32.A.Point), Tangent: cv(h32.A.Tangent)},
			B: f64hermit.PointTangent{Point: cv(h32.B.Point), Tangent: cv(h32.B.Tangent)},
		}
		b32 := h32.BoundingBox()
		b64 := h64.BoundingBox()
		for axis := 0; axis < 3; axis++ {
			relTol := 1e-5 * math.Max(math.Abs(b64.Min[axis]), math.Abs(b64.Max[axis]))
			if d := math.Abs(float64(b32.Min[axis]) - b64.Min[axis]); d > relTol {
				t.Fatalf("trial %d Min[%d]: f32 %.10g f64 %.10g diff %g tol %g", trial, axis, b32.Min[axis], b64.Min[axis], d, relTol)
			}
			if d := math.Abs(float64(b32.Max[axis]) - b64.Max[axis]); d > relTol {
				t.Fatalf("trial %d Max[%d]: f32 %.10g f64 %.10g diff %g tol %g", trial, axis, b32.Max[axis], b64.Max[axis], d, relTol)
			}
		}
	}
}
