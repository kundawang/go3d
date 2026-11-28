package vec3

import (
	"errors"
	"math"
)

// Triangle is a finite triangle in 3D space defined by its three vertices
// A, B and C.
type Triangle struct {
	A T
	B T
	C T
}

// TriangleRegion classifies the feature of a triangle that contains a closest
// point.
type TriangleRegion int

const (
	// RegionUnknown indicates that no closest point was computed.
	RegionUnknown TriangleRegion = iota

	// RegionInterior indicates a point in the interior of the triangle.
	RegionInterior
	// RegionEdgeAB indicates a point on the edge from A to B.
	RegionEdgeAB
	// RegionEdgeBC indicates a point on the edge from B to C.
	RegionEdgeBC
	// RegionEdgeCA indicates a point on the edge from C to A.
	RegionEdgeCA
	// RegionVertexA indicates the vertex A.
	RegionVertexA
	// RegionVertexB indicates the vertex B.
	RegionVertexB
	// RegionVertexC indicates the vertex C.
	RegionVertexC
)

// String returns the name of the region.
func (region TriangleRegion) String() string {
	switch region {
	case RegionInterior:
		return "interior"
	case RegionEdgeAB:
		return "edge AB"
	case RegionEdgeBC:
		return "edge BC"
	case RegionEdgeCA:
		return "edge CA"
	case RegionVertexA:
		return "vertex A"
	case RegionVertexB:
		return "vertex B"
	case RegionVertexC:
		return "vertex C"
	default:
		return "unknown"
	}
}

// ErrNonFiniteInput is returned by Triangle.ClosestPoint when the query point
// or a triangle vertex is NaN or infinite, or when a difference between two of
// them overflows.
var ErrNonFiniteInput = errors.New("vec3: closest point query requires finite input")

const (
	// machineEpsilon is the difference between 1 and the smallest float32
	// greater than 1.
	machineEpsilon = 1.1920928955078125e-07

	// triangleWeightEpsilon is the tolerated distance of a barycentric weight
	// from zero. Smaller weights are snapped to zero and the result is
	// reclassified. It is scaled with the machine epsilon so the same rule
	// applies to the float32 and float64 variants.
	triangleWeightEpsilon float32 = 16 * machineEpsilon

	// triangleDegenerateSinSqr is the relative threshold on the squared sine of
	// the angle between two edges. The determinant used to project into the
	// interior has a rounding noise of about machineEpsilon relative to the
	// product of the squared edge lengths, so a smaller squared sine carries no
	// information about the triangle's area. Below this threshold the interior
	// is ignored and the closest point is taken from the vertices and edges,
	// which form the triangle's convex hull.
	triangleDegenerateSinSqr float32 = 4 * machineEpsilon
)

// ClosestPointResult is the result of Triangle.ClosestPoint.
type ClosestPointResult struct {
	// Point is the closest point on the triangle to the query point.
	Point T

	// DistanceSqr is the squared distance between the query point and Point.
	// It is +Inf if the true squared distance is not representable.
	DistanceSqr float32

	// Weights holds the barycentric weights of the vertices A, B and C, in
	// that order. The weights are non-negative and sum to one up to floating
	// point precision.
	Weights T

	// Region tells if Point lies in the interior, on an edge or at a vertex.
	Region TriangleRegion
}

// ClosestPoint returns the point on the finite triangle closest to p.
//
// The query supports triangles with any orientation in 3D and points inside or
// outside the triangle. The result contains the closest point, its squared
// distance to p, the barycentric weights of the vertices A, B and C (in that
// order) and the region that contains the closest point.
//
// Weights are non-negative and sum to one up to floating point precision.
// Weights within triangleWeightEpsilon of zero are snapped to zero, the
// remaining weights are renormalized and the region is derived from the
// resulting zeros. This keeps the classification stable for points that lie on
// or very close to an edge or a vertex.
//
// Degenerate triangles are handled without an error:
//
//   - If all three vertices coincide, the closest point is that vertex with the
//     weights {1, 0, 0} and RegionVertexA.
//   - If the triangle has zero area (collinear vertices or two coincident
//     vertices), the interior is ignored and the closest point is the closest
//     point on the vertices and edges, which form the triangle's convex hull.
//
// Points are evaluated in the order A, B, C, AB, BC, CA, interior, and the
// first candidate with the smallest distance wins. The order makes the result
// deterministic for coincident vertices and for points that are equidistant
// from several features.
//
// Extremely small and large triangles are supported by scaling the three
// difference vectors to a common bounded frame before computing dot products
// and barycentric coordinates. This is a similarity transform, so it does not
// change any distance ratio. If the ratio between the query distance and the
// triangle size exceeds the float range, the result is still finite and stable,
// but it reflects no more than the leading direction.
//
// ErrNonFiniteInput is returned if the query point or a vertex is NaN or
// infinite, or if a difference between two of them overflows. In that case the
// result is the zero value with RegionUnknown.
func (tri *Triangle) ClosestPoint(p *T) (ClosestPointResult, error) {
	if !finiteVec(p) || !finiteVec(&tri.A) || !finiteVec(&tri.B) || !finiteVec(&tri.C) {
		return ClosestPointResult{Region: RegionUnknown}, ErrNonFiniteInput
	}

	a := tri.A
	ab := Sub(&tri.B, &a)
	ac := Sub(&tri.C, &a)
	ap := Sub(p, &a)
	if !finiteVec(&ab) || !finiteVec(&ac) || !finiteVec(&ap) {
		return ClosestPointResult{Region: RegionUnknown}, ErrNonFiniteInput
	}

	// Work in a similarity frame that scales all three difference vectors by
	// the same factor. This keeps dot products and barycentric coordinates in
	// a bounded range for extremely small or large triangles without changing
	// any distance ratio. Componentwise division preserves relative precision.
	scale := maxAbsComponent(&ab, &ac, &ap)
	if scale == 0 {
		// The vertices and the query point all coincide.
		return ClosestPointResult{
			Point:       a,
			DistanceSqr: 0,
			Weights:     T{1, 0, 0},
			Region:      RegionVertexA,
		}, nil
	}
	nab := divideVec(&ab, scale)
	nac := divideVec(&ac, scale)
	nap := divideVec(&ap, scale)

	// Candidates are ranked by their squared distance in the scaled frame,
	// where every component is at most one and no overflow can occur. The
	// returned point and distance are computed from the original coordinates.
	bestDist := float32(math.Inf(1))
	bestWeights := T{}
	haveBest := false

	// consider evaluates the candidate described by w. The point is computed
	// from the original vertices, so a vertex or edge is reproduced exactly.
	// The weights are normalized first so that the scaled-frame distance is
	// the true affine distance and the ranking is unaffected by rounding of
	// expressions such as 1-t+t.
	consider := func(w T) {
		sum := w[0] + w[1] + w[2]
		if sum > 0 && sum != 1 {
			w[0] /= sum
			w[1] /= sum
			w[2] /= sum
		}
		qn := T{
			w[1]*nab[0] + w[2]*nac[0],
			w[1]*nab[1] + w[2]*nac[1],
			w[1]*nab[2] + w[2]*nac[2],
		}
		dn := Sub(&nap, &qn)
		dist := dn.LengthSqr()
		if !haveBest || dist < bestDist {
			bestDist = dist
			bestWeights = w
			haveBest = true
		}
	}

	// Vertices, in deterministic order for coincident candidates.
	consider(T{1, 0, 0})
	consider(T{0, 1, 0})
	consider(T{0, 0, 1})

	// Edge AB.
	{
		d := Dot(&nab, &nab)
		t := float32(0)
		if d > 0 {
			t = clamp01(Dot(&nap, &nab) / d)
		}
		consider(T{1 - t, t, 0})
	}

	// Edge BC.
	{
		nbc := Sub(&nac, &nab)
		napb := Sub(&nap, &nab)
		d := Dot(&nbc, &nbc)
		t := float32(0)
		if d > 0 {
			t = clamp01(Dot(&napb, &nbc) / d)
		}
		consider(T{0, 1 - t, t})
	}

	// Edge CA, parameterized from C to A.
	{
		napc := Sub(&nap, &nac)
		d := Dot(&nac, &nac)
		t := float32(0)
		if d > 0 {
			t = clamp01(-Dot(&napc, &nac) / d)
		}
		consider(T{t, 0, 1 - t})
	}

	// Interior projection. It is only valid if the projection falls inside the
	// triangle. A zero or tiny cross product means that the triangle is
	// (numerically) degenerate; the vertices and edges above then already form
	// its convex hull.
	{
		d00 := Dot(&nab, &nab)
		d11 := Dot(&nac, &nac)
		d01 := Dot(&nab, &nac)
		denom := d00*d11 - d01*d01
		if denom > triangleDegenerateSinSqr*d00*d11 {
			d20 := Dot(&nap, &nab)
			d21 := Dot(&nap, &nac)
			v := (d11*d20 - d01*d21) / denom
			w := (d00*d21 - d01*d20) / denom
			u := 1 - v - w
			// The comparisons also reject NaN coordinates.
			if u >= 0 && v >= 0 && w >= 0 {
				consider(T{u, v, w})
			}
		}
	}

	if !haveBest {
		// Unreachable: the three vertices are always candidates.
		return ClosestPointResult{Region: RegionUnknown}, ErrNonFiniteInput
	}

	weights, region := canonicalWeights(bestWeights)
	q := weightedPoint(&tri.A, &tri.B, &tri.C, weights)
	d := Sub(p, &q)
	return ClosestPointResult{
		Point:       q,
		DistanceSqr: d.LengthSqr(),
		Weights:     weights,
		Region:      region,
	}, nil
}

// finiteVec reports whether all components of v are finite.
func finiteVec(v *T) bool {
	for _, x := range v {
		f := float64(x)
		if math.IsNaN(f) || math.IsInf(f, 0) {
			return false
		}
	}
	return true
}

// maxAbsComponent returns the largest absolute component of a, b and c.
func maxAbsComponent(a, b, c *T) float32 {
	m := abs32(a[0])
	for _, x := range [8]float32{a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]} {
		if abs := abs32(x); abs > m {
			m = abs
		}
	}
	return m
}

// abs32 returns the absolute value of x.
func abs32(x float32) float32 {
	if x < 0 {
		return -x
	}
	return x
}

// divideVec returns v divided by s.
func divideVec(v *T, s float32) T {
	return T{v[0] / s, v[1] / s, v[2] / s}
}

// clamp01 clamps t to [0, 1]. NaN maps to 0, the deterministic fallback for
// derived values that overflowed.
func clamp01(t float32) float32 {
	switch {
	case t > 1:
		return 1
	case t > 0:
		return t
	default:
		return 0
	}
}

// weightedPoint returns the affine combination w[0]*a + w[1]*b + w[2]*c.
func weightedPoint(a, b, c *T, w T) T {
	return T{
		w[0]*a[0] + w[1]*b[0] + w[2]*c[0],
		w[0]*a[1] + w[1]*b[1] + w[2]*c[1],
		w[0]*a[2] + w[1]*b[2] + w[2]*c[2],
	}
}

// canonicalWeights snaps, renormalizes and clamps barycentric weights and
// derives the region from the resulting zeros.
func canonicalWeights(w T) (T, TriangleRegion) {
	for i := range w {
		// This also turns NaN into zero.
		if !(w[i] > triangleWeightEpsilon) {
			w[i] = 0
		}
	}

	zeros := 0
	zeroIndex := 0
	for i, x := range w {
		if x == 0 {
			zeros++
			zeroIndex = i
		}
	}

	if zeros >= 2 {
		// A vertex. If rounding left more than one non-zero weight, keep the
		// largest one; ties keep the first vertex in A, B, C order.
		best := 0
		for i := 1; i < 3; i++ {
			if w[i] > w[best] {
				best = i
			}
		}
		var out T
		out[best] = 1
		return out, vertexRegion(best)
	}

	sum := w[0] + w[1] + w[2]
	if !(sum > 0) {
		// Cannot happen for a valid candidate; keep the classification total.
		return T{1, 0, 0}, RegionVertexA
	}
	for i := range w {
		w[i] /= sum
	}

	if zeros == 1 {
		return w, edgeRegion(zeroIndex)
	}
	return w, RegionInterior
}

// vertexRegion returns the region of the vertex with the given weight index.
func vertexRegion(index int) TriangleRegion {
	switch index {
	case 0:
		return RegionVertexA
	case 1:
		return RegionVertexB
	default:
		return RegionVertexC
	}
}

// edgeRegion returns the edge opposite to the vertex with the given weight
// index, which is the edge that remains when that weight is zero.
func edgeRegion(zeroIndex int) TriangleRegion {
	switch zeroIndex {
	case 0:
		return RegionEdgeBC
	case 1:
		return RegionEdgeCA
	default:
		return RegionEdgeAB
	}
}
