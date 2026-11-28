package vec3

import (
	"math"
	"testing"
)

// rightTriangle is a unit right triangle in the z=0 plane.
var rightTriangle = Triangle{
	A: T{0, 0, 0},
	B: T{1, 0, 0},
	C: T{0, 1, 0},
}

func checkClosest(t *testing.T, name string, tri Triangle, p T, wantPoint T, wantDistSqr float64, wantWeights T, wantRegion TriangleRegion) {
	t.Helper()

	got, err := tri.ClosestPoint(&p)
	if err != nil {
		t.Fatalf("%s: unexpected error: %v", name, err)
	}

	if !got.Point.PracticallyEquals(&wantPoint, 1e-12) {
		t.Errorf("%s: point = %v, want %v", name, got.Point, wantPoint)
	}
	if math.Abs(got.DistanceSqr-wantDistSqr) > 1e-12*math.Max(1, math.Abs(wantDistSqr)) {
		t.Errorf("%s: distance sqr = %v, want %v", name, got.DistanceSqr, wantDistSqr)
	}
	if !got.Weights.PracticallyEquals(&wantWeights, 1e-12) {
		t.Errorf("%s: weights = %v, want %v", name, got.Weights, wantWeights)
	}
	if got.Region != wantRegion {
		t.Errorf("%s: region = %v, want %v", name, got.Region, wantRegion)
	}

	// Invariants that must hold for every result.
	sum := got.Weights[0] + got.Weights[1] + got.Weights[2]
	if math.Abs(sum-1) > 1e-12 {
		t.Errorf("%s: weights sum to %v, want 1", name, sum)
	}
	for i, w := range got.Weights {
		if w < 0 {
			t.Errorf("%s: weight %d = %v is negative", name, i, w)
		}
	}
	reconstructed := weightedPoint(&tri.A, &tri.B, &tri.C, got.Weights)
	if !reconstructed.PracticallyEquals(&got.Point, 1e-12) {
		t.Errorf("%s: weights reconstruct %v, point is %v", name, reconstructed, got.Point)
	}
}

func TestTriangleClosestPointRightTriangle(t *testing.T) {
	cases := []struct {
		name    string
		p       T
		point   T
		distSqr float64
		weights T
		region  TriangleRegion
	}{
		{"above interior", T{0.25, 0.25, 2}, T{0.25, 0.25, 0}, 4, T{0.5, 0.25, 0.25}, RegionInterior},
		{"below interior", T{0.25, 0.25, -3}, T{0.25, 0.25, 0}, 9, T{0.5, 0.25, 0.25}, RegionInterior},
		{"on interior", T{0.25, 0.25, 0}, T{0.25, 0.25, 0}, 0, T{0.5, 0.25, 0.25}, RegionInterior},
		{"vertex A", T{-3, -4, 0}, T{0, 0, 0}, 25, T{1, 0, 0}, RegionVertexA},
		{"vertex B", T{5, 0, 0}, T{1, 0, 0}, 16, T{0, 1, 0}, RegionVertexB},
		{"vertex C", T{0, 5, 0}, T{0, 1, 0}, 16, T{0, 0, 1}, RegionVertexC},
		{"edge AB", T{0.5, 0, 0}, T{0.5, 0, 0}, 0, T{0.5, 0.5, 0}, RegionEdgeAB},
		{"edge AB outside", T{0.5, -3, 0}, T{0.5, 0, 0}, 9, T{0.5, 0.5, 0}, RegionEdgeAB},
		{"edge BC", T{0.5, 0.5, 0}, T{0.5, 0.5, 0}, 0, T{0, 0.5, 0.5}, RegionEdgeBC},
		{"edge BC outside", T{2, 2, 0}, T{0.5, 0.5, 0}, 4.5, T{0, 0.5, 0.5}, RegionEdgeBC},
		{"edge CA", T{0, 0.5, 0}, T{0, 0.5, 0}, 0, T{0.5, 0, 0.5}, RegionEdgeCA},
		{"outside near B corner", T{2, -1, 0}, T{1, 0, 0}, 2, T{0, 1, 0}, RegionVertexB},
		{"outside beyond C corner", T{-1, 2, 0}, T{0, 1, 0}, 2, T{0, 0, 1}, RegionVertexC},
		{"outside beyond A corner", T{-1, -1, 0}, T{0, 0, 0}, 2, T{1, 0, 0}, RegionVertexA},
	}

	for _, c := range cases {
		checkClosest(t, c.name, rightTriangle, c.p, c.point, c.distSqr, c.weights, c.region)
	}
}

func TestTriangleClosestPointAnyOrientation(t *testing.T) {
	// A triangle not aligned with any coordinate plane.
	tri := Triangle{
		A: T{0, 0, 0},
		B: T{1, 1, 0},
		C: T{-1, 1, 1},
	}
	// The centroid projected 2 away along the triangle normal (1, -1, 2).
	// The squared length of the normal is 6.
	centroid := weightedPoint(&tri.A, &tri.B, &tri.C, T{1.0 / 3, 1.0 / 3, 1.0 / 3})
	ab := Sub(&tri.B, &tri.A)
	ac := Sub(&tri.C, &tri.A)
	normal := Cross(&ab, &ac)
	p := centroid
	scaled := normal.Scaled(2)
	p.Add(&scaled)

	checkClosest(t, "along normal", tri, p, centroid, 24, T{1.0 / 3, 1.0 / 3, 1.0 / 3}, RegionInterior)

	// A point on the far side of the plane behaves the same way.
	p = centroid
	scaled = normal.Scaled(3)
	p.Sub(&scaled)
	checkClosest(t, "opposite normal", tri, p, centroid, 54, T{1.0 / 3, 1.0 / 3, 1.0 / 3}, RegionInterior)

	// Brute force check that no point sampled on the triangle is closer.
	p = T{3, -2, 4}
	got, err := tri.ClosestPoint(&p)
	if err != nil {
		t.Fatal(err)
	}
	const steps = 50
	for i := 0; i <= steps; i++ {
		for j := 0; i+j <= steps; j++ {
			v := float64(i) / steps
			w := float64(j) / steps
			sample := weightedPoint(&tri.A, &tri.B, &tri.C, T{1 - v - w, v, w})
			diff := Sub(&p, &sample)
			if d := diff.LengthSqr(); d < got.DistanceSqr-1e-9 {
				t.Fatalf("sample %v is closer (%v) than reported closest point (%v)", sample, d, got.DistanceSqr)
			}
		}
	}
}

func TestTriangleClosestPointDegenerate(t *testing.T) {
	t.Run("all vertices coincident", func(t *testing.T) {
		v := T{1, 2, 3}
		tri := Triangle{A: v, B: v, C: v}
		checkClosest(t, "coincident", tri, T{1, 2, 4}, v, 1, T{1, 0, 0}, RegionVertexA)
	})

	t.Run("collinear horizontal", func(t *testing.T) {
		tri := Triangle{A: T{0, 0, 0}, B: T{1, 0, 0}, C: T{2, 0, 0}}
		// The closest point is on the hull segment from x=0 to x=2.
		checkClosest(t, "collinear", tri, T{1.5, 1, 0}, T{1.5, 0, 0}, 1, T{0, 0.5, 0.5}, RegionEdgeBC)
	})

	t.Run("collinear vertical", func(t *testing.T) {
		tri := Triangle{A: T{0, 0, 0}, B: T{0, 0, 1}, C: T{0, 0, 2}}
		checkClosest(t, "collinear", tri, T{1, 0, 1.5}, T{0, 0, 1.5}, 1, T{0, 0.5, 0.5}, RegionEdgeBC)
	})

	t.Run("duplicate A and B", func(t *testing.T) {
		tri := Triangle{A: T{0, 0, 0}, B: T{0, 0, 0}, C: T{0, 1, 0}}
		checkClosest(t, "duplicate", tri, T{1, 0.5, 0}, T{0, 0.5, 0}, 1, T{0, 0.5, 0.5}, RegionEdgeBC)
	})

	t.Run("duplicate B and C", func(t *testing.T) {
		tri := Triangle{A: T{0, 0, 0}, B: T{0, 1, 0}, C: T{0, 1, 0}}
		// B and C coincide. The deterministic candidate order classifies the
		// hit as vertex B instead of the zero-length edge BC.
		checkClosest(t, "duplicate", tri, T{0, 3, 0}, T{0, 1, 0}, 4, T{0, 1, 0}, RegionVertexB)
	})

	t.Run("nearly collinear is treated as degenerate", func(t *testing.T) {
		tri := Triangle{A: T{0, 0, 0}, B: T{1, 0, 0}, C: T{2, 1e-15, 0}}
		got, err := tri.ClosestPoint(&T{1.5, 1, 0})
		if err != nil {
			t.Fatal(err)
		}
		if got.DistanceSqr > 1+1e-6 {
			t.Errorf("distance sqr = %v, want about 1", got.DistanceSqr)
		}
		// The projection point lies on the nearly flat triangle, but the
		// interior is not resolvable, so the result is on the hull.
		if got.Region == RegionInterior {
			t.Errorf("region = interior, want a hull feature")
		}
	})
}

func TestTriangleClosestPointExtremeScale(t *testing.T) {
	t.Run("tiny", func(t *testing.T) {
		const s = 1e-20
		tri := Triangle{A: T{0, 0, 0}, B: T{s, 0, 0}, C: T{0, s, 0}}
		p := T{0.25 * s, 0.25 * s, s}
		got, err := tri.ClosestPoint(&p)
		if err != nil {
			t.Fatal(err)
		}
		want := T{0.25 * s, 0.25 * s, 0}
		if !got.Point.PracticallyEquals(&want, 1e-30) {
			t.Errorf("point = %v, want %v", got.Point, want)
		}
		if !got.Weights.PracticallyEquals(&T{0.5, 0.25, 0.25}, 1e-12) {
			t.Errorf("weights = %v", got.Weights)
		}
		if got.Region != RegionInterior {
			t.Errorf("region = %v, want interior", got.Region)
		}
		if rel := math.Abs(got.DistanceSqr-s*s) / (s * s); rel > 1e-9 {
			t.Errorf("distance sqr = %v, want about %v", got.DistanceSqr, s*s)
		}
	})

	t.Run("large", func(t *testing.T) {
		const s = 1e300
		tri := Triangle{A: T{0, 0, 0}, B: T{s, 0, 0}, C: T{0, s, 0}}
		p := T{0.25 * s, 0.25 * s, s}
		got, err := tri.ClosestPoint(&p)
		if err != nil {
			t.Fatal(err)
		}
		// The squared distance is not representable, it must saturate at +Inf
		// rather than become NaN.
		if !math.IsInf(float64(got.DistanceSqr), 1) {
			t.Errorf("distance sqr = %v, want +Inf", got.DistanceSqr)
		}
		if !got.Weights.PracticallyEquals(&T{0.5, 0.25, 0.25}, 1e-12) {
			t.Errorf("weights = %v", got.Weights)
		}
		if got.Region != RegionInterior {
			t.Errorf("region = %v, want interior", got.Region)
		}
	})

	t.Run("far point tiny triangle", func(t *testing.T) {
		// A finite point whose ratio to the triangle size exceeds the float
		// range cannot be resolved; the query must stay finite and stable.
		tri := Triangle{A: T{0, 0, 0}, B: T{1e-300, 0, 0}, C: T{0, 1e-300, 0}}
		p := T{1e300, 1e300, 1e300}
		got, err := tri.ClosestPoint(&p)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if math.IsNaN(float64(got.DistanceSqr)) {
			t.Errorf("distance sqr is NaN")
		}
		sum := got.Weights[0] + got.Weights[1] + got.Weights[2]
		if math.Abs(sum-1) > 1e-12 {
			t.Errorf("weights %v do not sum to 1", got.Weights)
		}
		if got.Region == RegionUnknown {
			t.Errorf("region is unknown")
		}
	})
}

func TestTriangleClosestPointNonFinite(t *testing.T) {
	nan := math.NaN()
	inf := math.Inf(1)
	ninf := math.Inf(-1)

	t.Run("query point", func(t *testing.T) {
		for _, p := range []T{{nan, 0, 0}, {0, inf, 0}, {0, 0, ninf}} {
			got, err := rightTriangle.ClosestPoint(&p)
			if err != ErrNonFiniteInput {
				t.Errorf("point %v: err = %v, want %v", p, err, ErrNonFiniteInput)
			}
			if got.Region != RegionUnknown {
				t.Errorf("point %v: region = %v, want unknown", p, got.Region)
			}
		}
	})

	t.Run("vertices", func(t *testing.T) {
		p := T{0.25, 0.25, 1}
		bad := []Triangle{
			{A: T{nan, 0, 0}, B: T{1, 0, 0}, C: T{0, 1, 0}},
			{A: T{0, 0, 0}, B: T{inf, 0, 0}, C: T{0, 1, 0}},
			{A: T{0, 0, 0}, B: T{1, 0, 0}, C: T{0, 0, ninf}},
		}
		for i, tri := range bad {
			if _, err := tri.ClosestPoint(&p); err != ErrNonFiniteInput {
				t.Errorf("triangle %d: err = %v, want %v", i, err, ErrNonFiniteInput)
			}
		}
	})

	t.Run("overflowing difference", func(t *testing.T) {
		tri := Triangle{A: T{-math.MaxFloat64, 0, 0}, B: T{math.MaxFloat64, 0, 0}, C: T{0, 1, 0}}
		if _, err := tri.ClosestPoint(&T{0, 0, 1}); err != ErrNonFiniteInput {
			t.Errorf("err = %v, want %v", err, ErrNonFiniteInput)
		}
	})
}

func TestTriangleWeightsInterpolateUV(t *testing.T) {
	// Per-vertex UV coordinates, as used by a level editor decal.
	uvA := [2]float64{0, 0}
	uvB := [2]float64{1, 0}
	uvC := [2]float64{0, 1}

	got, err := rightTriangle.ClosestPoint(&T{0.25, 0.25, 3})
	if err != nil {
		t.Fatal(err)
	}

	uv := [2]float64{
		got.Weights[0]*uvA[0] + got.Weights[1]*uvB[0] + got.Weights[2]*uvC[0],
		got.Weights[0]*uvA[1] + got.Weights[1]*uvB[1] + got.Weights[2]*uvC[1],
	}
	want := [2]float64{0.25, 0.25}
	if math.Abs(uv[0]-want[0]) > 1e-12 || math.Abs(uv[1]-want[1]) > 1e-12 {
		t.Errorf("interpolated uv = %v, want %v", uv, want)
	}
}

func TestTriangleRegionString(t *testing.T) {
	cases := map[TriangleRegion]string{
		RegionUnknown:  "unknown",
		RegionInterior: "interior",
		RegionEdgeAB:   "edge AB",
		RegionEdgeBC:   "edge BC",
		RegionEdgeCA:   "edge CA",
		RegionVertexA:  "vertex A",
		RegionVertexB:  "vertex B",
		RegionVertexC:  "vertex C",
	}
	for region, want := range cases {
		if got := region.String(); got != want {
			t.Errorf("%d: String() = %q, want %q", int(region), got, want)
		}
	}
}
