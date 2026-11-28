package vec3_test

import (
	"math"
	"math/rand"
	"testing"

	"github.com/ungerik/go3d/float64/vec3"
	vec3f32 "github.com/ungerik/go3d/vec3"
)

// TestTriangleClosestPointCrossPrecision checks that the float32 and float64
// implementations agree on well-conditioned inputs. Borderline hits, where the
// closest point lies within float32 resolution of an edge or vertex, are
// excluded because no float32 implementation can resolve them.
func TestTriangleClosestPointCrossPrecision(t *testing.T) {
	rng := rand.New(rand.NewSource(1))

	for i := 0; i < 2000; i++ {
		// Random, reasonably conditioned triangles and points.
		vals := make([]float64, 12)
		for j := range vals {
			vals[j] = (rng.Float64()*2 - 1) * 10
		}
		tri64 := vec3.Triangle{
			A: vec3.T{vals[0], vals[1], vals[2]},
			B: vec3.T{vals[3], vals[4], vals[5]},
			C: vec3.T{vals[6], vals[7], vals[8]},
		}
		tri32 := vec3f32.Triangle{
			A: vec3f32.T{float32(vals[0]), float32(vals[1]), float32(vals[2])},
			B: vec3f32.T{float32(vals[3]), float32(vals[4]), float32(vals[5])},
			C: vec3f32.T{float32(vals[6]), float32(vals[7]), float32(vals[8])},
		}

		// Build the query point from barycentric weights. Points close to a
		// feature are skipped so only resolvable classifications are compared.
		u := 0.05 + 0.9*rng.Float64()
		v := 0.05 + (0.9-u)*rng.Float64()
		w := 1 - u - v
		if w < 0.05 {
			continue
		}
		raw := vec3.T{
			u*tri64.A[0] + v*tri64.B[0] + w*tri64.C[0],
			u*tri64.A[1] + v*tri64.B[1] + w*tri64.C[1],
			u*tri64.A[2] + v*tri64.B[2] + w*tri64.C[2],
		}
		p64 := raw
		p64[2] += rng.Float64()*5 + 0.5
		p32 := vec3f32.T{float32(p64[0]), float32(p64[1]), float32(p64[2])}

		got64, err64 := tri64.ClosestPoint(&p64)
		got32, err32 := tri32.ClosestPoint(&p32)
		if (err64 == nil) != (err32 == nil) {
			t.Fatalf("case %d: error mismatch %v / %v", i, err64, err32)
		}
		if err64 != nil {
			continue
		}

		if got32.Region != vec3f32.TriangleRegion(got64.Region) {
			t.Errorf("case %d: region %v (float32) vs %v (float64)", i, got32.Region, got64.Region)
			continue
		}
		for k := 0; k < 3; k++ {
			if math.Abs(float64(got32.Weights[k])-got64.Weights[k]) > 1e-4 {
				t.Errorf("case %d: weight %d = %v (float32) vs %v (float64)", i, k, got32.Weights[k], got64.Weights[k])
			}
		}
		d32 := math.Sqrt(float64(got32.DistanceSqr))
		d64 := math.Sqrt(got64.DistanceSqr)
		if math.Abs(d32-d64) > 1e-4*math.Max(1, d64) {
			t.Errorf("case %d: distance %v (float32) vs %v (float64)", i, d32, d64)
		}
	}
}

// TestTriangleClosestPointCrossPrecisionFeatures checks exact feature hits:
// vertices and edge midpoints must be classified identically by both
// precisions.
func TestTriangleClosestPointCrossPrecisionFeatures(t *testing.T) {
	tri64 := vec3.Triangle{
		A: vec3.T{0.125, -2.5, 3},
		B: vec3.T{4, 1.75, -0.5},
		C: vec3.T{-1.5, 6, 2.25},
	}
	tri32 := vec3f32.Triangle{
		A: vec3f32.T{0.125, -2.5, 3},
		B: vec3f32.T{4, 1.75, -0.5},
		C: vec3f32.T{-1.5, 6, 2.25},
	}

	mid := func(a, b vec3.T) vec3.T {
		return vec3.T{(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2}
	}

	// Offset along the plane normal so the orthogonal projection of the lifted
	// point is exactly the feature on the triangle.
	ab := vec3.T{tri64.B[0] - tri64.A[0], tri64.B[1] - tri64.A[1], tri64.B[2] - tri64.A[2]}
	ac := vec3.T{tri64.C[0] - tri64.A[0], tri64.C[1] - tri64.A[1], tri64.C[2] - tri64.A[2]}
	normal := vec3.Cross(&ab, &ac)
	normal.Normalize()
	offset := normal.Scaled(3)

	cases := []struct {
		name   string
		p      vec3.T
		region vec3.TriangleRegion
	}{
		{"vertex A", tri64.A, vec3.RegionVertexA},
		{"vertex B", tri64.B, vec3.RegionVertexB},
		{"vertex C", tri64.C, vec3.RegionVertexC},
		{"edge AB", mid(tri64.A, tri64.B), vec3.RegionEdgeAB},
		{"edge BC", mid(tri64.B, tri64.C), vec3.RegionEdgeBC},
		{"edge CA", mid(tri64.C, tri64.A), vec3.RegionEdgeCA},
	}

	for _, c := range cases {
		p64 := c.p
		p64.Add(&offset)
		p32 := vec3f32.T{float32(p64[0]), float32(p64[1]), float32(p64[2])}

		got64, err64 := tri64.ClosestPoint(&p64)
		got32, err32 := tri32.ClosestPoint(&p32)
		if err64 != nil || err32 != nil {
			t.Fatalf("%s: unexpected errors %v / %v", c.name, err64, err32)
		}
		if got64.Region != c.region {
			t.Errorf("%s: float64 region = %v, want %v", c.name, got64.Region, c.region)
		}
		if got32.Region != vec3f32.TriangleRegion(c.region) {
			t.Errorf("%s: float32 region = %v, want %v", c.name, got32.Region, c.region)
		}
	}
}
