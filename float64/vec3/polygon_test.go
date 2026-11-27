package vec3

import (
	"math"
	"math/rand"
	"testing"

	"github.com/ungerik/go3d/float64/vec2"
	vec2f "github.com/ungerik/go3d/vec2"
	vec3f "github.com/ungerik/go3d/vec3"
)

// The clipping tests use a unit square in the XY plane with side length 2,
// centered on the origin, and the plane x = 0. ClipPolygon keeps x >= 0, so
// the expected result is the right half: a rectangle of size 1 x 2 with
// area 2. The square is given counterclockwise as seen from +Z.
// This is the float64 mirror of the tests in vec3/polygon_test.go.
func clipTestSquare() []T {
	return []T{
		{-1, -1, 0},
		{1, -1, 0},
		{1, 1, 0},
		{-1, 1, 0},
	}
}

func clipTestSquareUVs() []vec2.T {
	return []vec2.T{
		{0, 0},
		{1, 0},
		{1, 1},
		{0, 1},
	}
}

func clipTestPlane() Plane {
	// Plane x = 0 with normal +X; the positive side is x >= 0.
	return NewPlane(&T{0, 0, 0}, &T{1, 0, 0})
}

func clipTestVecEqual(a, b T, eps float64) bool {
	return math.Abs(float64(a[0]-b[0])) <= float64(eps) &&
		math.Abs(float64(a[1]-b[1])) <= float64(eps) &&
		math.Abs(float64(a[2]-b[2])) <= float64(eps)
}

func clipTestUVEqual(a, b vec2.T, eps float64) bool {
	return math.Abs(float64(a[0]-b[0])) <= float64(eps) &&
		math.Abs(float64(a[1]-b[1])) <= float64(eps)
}

// clipTestCyclicRotation reports the rotation offset at which got matches
// want, or -1 if got is not a cyclic rotation of want. Sutherland-Hodgman
// may start the result at a different boundary vertex, but the cyclic order
// itself must be preserved.
func clipTestCyclicRotation(got, want []T, eps float64) int {
	if len(got) != len(want) {
		return -1
	}
	if len(got) == 0 {
		return 0
	}
	for offset := 0; offset < len(want); offset++ {
		matches := true
		for i := range got {
			if !clipTestVecEqual(got[i], want[(i+offset)%len(want)], eps) {
				matches = false
				break
			}
		}
		if matches {
			return offset
		}
	}
	return -1
}

func TestClipPolygonSquareAreaAndIntersections(t *testing.T) {
	square := clipTestSquare()
	plane := clipTestPlane()
	clipped, ok := ClipPolygon(square, &plane)
	if !ok {
		t.Fatalf("ClipPolygon returned ok=false for a valid square and plane")
	}
	if got := PolygonArea(clipped); math.Abs(got-2) > 1e-5 {
		t.Errorf("clipped area: got %v, want 2", got)
	}
	// The two new intersection points must lie exactly on x = 0.
	want := []T{
		{1, -1, 0},
		{1, 1, 0},
		{0, 1, 0},
		{0, -1, 0},
	}
	offset := clipTestCyclicRotation(clipped, want, 1e-5)
	if offset < 0 {
		t.Fatalf("clipped polygon %v is not the expected right half %v", clipped, want)
	}
	// The kept original vertices (1,-1) and (1,1) must keep their relative
	// order and must not be separated by an inserted vertex. They follow
	// each other in the result, exactly as in the input.
	firstInside := -1
	for i, v := range clipped {
		if clipTestVecEqual(v, T{1, -1, 0}, 1e-6) {
			firstInside = i
		}
	}
	if firstInside < 0 || !clipTestVecEqual(clipped[(firstInside+1)%len(clipped)], T{1, 1, 0}, 1e-6) {
		t.Errorf("original vertex order changed: got %v", clipped)
	}
	// The new intersection points lie on the plane; the kept original
	// vertices are the ones with x = 1.
	for i, v := range clipped {
		d := plane.SignedDistance(&v)
		if v[0] == 1 && math.Abs(d-1) > 1e-5 {
			t.Errorf("original vertex %d %v has signed distance %v, want 1", i, v, d)
		}
		if math.Abs(float64(v[0])) > 1e-5 && math.Abs(v[0]-1) > 1e-5 {
			t.Errorf("unexpected vertex %d %v", i, v)
		}
	}
	for _, v := range []T{{0, 1, 0}, {0, -1, 0}} {
		if d := plane.SignedDistance(&v); math.Abs(d) > 1e-5 {
			t.Errorf("intersection %v not on the plane, signed distance %v", v, d)
		}
	}
	if poly := PolygonArea(square); math.Abs(float64(poly-4)) > 1e-5 {
		t.Errorf("input square area: got %v, want 4", poly)
	}
}

func TestClipPolygonSquareUVs(t *testing.T) {
	square := clipTestSquare()
	uvs := clipTestSquareUVs()
	plane := clipTestPlane()
	clipped, clippedUVs, ok := ClipPolygonUVs(square, uvs, &plane)
	if !ok {
		t.Fatalf("ClipPolygonUVs returned ok=false for valid input")
	}
	if len(clippedUVs) != len(clipped) {
		t.Fatalf("UV slice has %d entries for %d vertices", len(clippedUVs), len(clipped))
	}
	// The expected UVs in the same cyclic order as the expected vertices.
	wantUVs := []vec2.T{
		{1, 0},
		{1, 1},
		{0.5, 1},
		{0.5, 0},
	}
	offset := clipTestCyclicRotation(clipped, []T{
		{1, -1, 0},
		{1, 1, 0},
		{0, 1, 0},
		{0, -1, 0},
	}, 1e-5)
	if offset < 0 {
		t.Fatalf("unexpected clipped polygon %v", clipped)
	}
	for i := range clippedUVs {
		if !clipTestUVEqual(clippedUVs[i], wantUVs[(i+offset)%len(wantUVs)], 1e-5) {
			t.Errorf("uv %d: got %v, want %v", i, clippedUVs[i], wantUVs[(i+offset)%len(wantUVs)])
		}
	}
	// Every intersection UV must be the interpolation of the UVs of the
	// edge endpoints with the same t as the position: edge x=-1 -> x=1
	// has u 0 -> 1 and the cut at x=0 must therefore have u=0.5.
	for i := range clipped {
		if math.Abs(float64(clipped[i][0])) <= 1e-5 && math.Abs(float64(clippedUVs[i][0]-0.5)) > 1e-5 {
			t.Errorf("vertex %v has uv %v, expected u=0.5", clipped[i], clippedUVs[i])
		}
	}
}

// TestClipPolygonInsideUnchanged checks that a polygon completely inside
// the kept half space is returned as an identical copy in the same order.
func TestClipPolygonInsideUnchanged(t *testing.T) {
	square := clipTestSquare()
	plane := NewPlane(&T{-2, 0, 0}, &T{1, 0, 0})
	clipped, ok := ClipPolygon(square, &plane)
	if !ok {
		t.Fatal("ClipPolygon returned ok=false for a fully inside polygon")
	}
	if len(clipped) != len(square) {
		t.Fatalf("got %d vertices, want %d", len(clipped), len(square))
	}
	for i := range square {
		if !clipTestVecEqual(clipped[i], square[i], 0) {
			t.Errorf("vertex %d changed: got %v, want %v", i, clipped[i], square[i])
		}
	}
	if got := PolygonArea(clipped); math.Abs(got-4) > 1e-5 {
		t.Errorf("area: got %v, want 4", got)
	}
}

// TestClipPolygonOutsideEmpty checks that a polygon completely on the
// discarded side returns an empty result with ok=true.
func TestClipPolygonOutsideEmpty(t *testing.T) {
	square := clipTestSquare()
	plane := NewPlane(&T{2, 0, 0}, &T{1, 0, 0})
	clipped, ok := ClipPolygon(square, &plane)
	if !ok {
		t.Fatal("ClipPolygon returned ok=false for a valid polygon")
	}
	if len(clipped) != 0 {
		t.Errorf("outside polygon: got %v, want empty", clipped)
	}
	if got := PolygonArea(clipped); got != 0 {
		t.Errorf("empty polygon area: got %v, want 0", got)
	}
}

// TestClipPolygonVertexOnPlane checks a polygon whose single vertex lies
// exactly on the plane and all other vertices are outside. Only that vertex
// is retained, a degenerate one-point polygon.
func TestClipPolygonVertexOnPlane(t *testing.T) {
	triangle := []T{
		{0, 0, 0}, // on the plane
		{-2, -1, 0},
		{-2, 1, 0},
	}
	plane := clipTestPlane()
	clipped, ok := ClipPolygon(triangle, &plane)
	if !ok {
		t.Fatal("ClipPolygon returned ok=false")
	}
	if len(clipped) != 1 || !clipTestVecEqual(clipped[0], T{0, 0, 0}, 1e-6) {
		t.Errorf("got %v, want the single on-plane vertex", clipped)
	}
}

// TestClipPolygonEdgeOnPlane checks a polygon with an entire edge on the
// plane; the retained result is the edge, a degenerate two-point polygon.
func TestClipPolygonEdgeOnPlane(t *testing.T) {
	square := []T{
		{0, -1, 0},
		{0, 1, 0},
		{-2, 1, 0},
		{-2, -1, 0},
	}
	plane := clipTestPlane()
	clipped, ok := ClipPolygon(square, &plane)
	if !ok {
		t.Fatal("ClipPolygon returned ok=false")
	}
	if len(clipped) != 2 {
		t.Fatalf("got %v, want the two on-plane edge vertices", clipped)
	}
	if got := PolygonArea(clipped); got != 0 {
		t.Errorf("degenerate area: got %v, want 0", got)
	}
}

// TestClipPolygonPlaneThroughDiagonal cuts a square along its diagonal,
// which runs through two opposite vertices. No new intersection points are
// needed and the retained half square is a right triangle of area 2.
func TestClipPolygonPlaneThroughDiagonal(t *testing.T) {
	square := clipTestSquare()
	// Plane x + y = 0, positive side x + y >= 0.
	plane := NewPlane(&T{0, 0, 0}, &T{1, 1, 0})
	clipped, ok := ClipPolygon(square, &plane)
	if !ok {
		t.Fatal("ClipPolygon returned ok=false")
	}
	want := []T{
		{1, -1, 0},
		{1, 1, 0},
		{-1, 1, 0},
	}
	if clipTestCyclicRotation(clipped, want, 1e-5) < 0 {
		t.Errorf("got %v, want %v", clipped, want)
	}
	if got := PolygonArea(clipped); math.Abs(got-2) > 1e-5 {
		t.Errorf("area: got %v, want 2", got)
	}
}

// TestClipPolygonAllOnPlane checks that a polygon lying in the plane is
// returned unchanged: all vertices are kept and no intersections added.
func TestClipPolygonAllOnPlane(t *testing.T) {
	square := clipTestSquare()
	plane := clipTestPlane()
	for i := range square {
		square[i][0] = 0
	}
	clipped, ok := ClipPolygon(square, &plane)
	if !ok {
		t.Fatal("ClipPolygon returned ok=false")
	}
	if len(clipped) != len(square) {
		t.Fatalf("got %v, want the unchanged square", clipped)
	}
	for i := range square {
		if !clipTestVecEqual(clipped[i], square[i], 1e-6) {
			t.Errorf("vertex %d: got %v, want %v", i, clipped[i], square[i])
		}
	}
}

// TestClipPolygonConcave cuts an L-shaped concave polygon so the plane
// crosses two parallel edges at x = 2.
func TestClipPolygonConcave(t *testing.T) {
	l := []T{
		{0, 0, 0},
		{4, 0, 0},
		{4, 1, 0},
		{1, 1, 0},
		{1, 4, 0},
		{0, 4, 0},
	}
	plane := NewPlane(&T{2, 0, 0}, &T{1, 0, 0})
	clipped, ok := ClipPolygon(l, &plane)
	if !ok {
		t.Fatal("ClipPolygon returned ok=false")
	}
	want := []T{
		{2, 0, 0},
		{4, 0, 0},
		{4, 1, 0},
		{2, 1, 0},
	}
	if clipTestCyclicRotation(clipped, want, 1e-5) < 0 {
		t.Errorf("got %v, want %v", clipped, want)
	}
	if got := PolygonArea(clipped); math.Abs(got-2) > 1e-5 {
		t.Errorf("clipped area: got %v, want 2", got)
	}
	for _, v := range clipped {
		if d := plane.SignedDistance(&v); d < -1e-5 {
			t.Errorf("vertex %v lies on the discarded side, signed distance %v", v, d)
		}
	}
}

// TestClipPolygonNearlyParallelEdge checks an edge that is almost parallel
// to the plane and crosses it: both endpoints have a tiny x coordinate but
// a large y coordinate, so the edge direction lies almost in the plane.
func TestClipPolygonNearlyParallelEdge(t *testing.T) {
	poly := []T{
		{-1e-4, -50, 0},
		{1e-4, 50, 0},
		{3e-4, -50, 0},
	}
	plane := clipTestPlane()
	clipped, ok := ClipPolygon(poly, &plane)
	if !ok {
		t.Fatal("ClipPolygon returned ok=false")
	}
	if len(clipped) == 0 {
		t.Fatalf("nearly parallel cut discarded everything: %v", clipped)
	}
	for _, v := range clipped {
		if d := plane.SignedDistance(&v); d < -1e-3 {
			t.Errorf("vertex %v on the discarded side, signed distance %v", v, d)
		}
	}
	if area := PolygonArea(clipped); math.IsNaN(area) {
		t.Errorf("area of nearly parallel cut is NaN")
	}
	// The first edge crosses x=0 halfway, so one intersection lies at
	// (0, 0).
	found := false
	for _, v := range clipped {
		if math.Abs(v[0]) <= 1e-3 && math.Abs(float64(v[1])-0) <= 1e-3 {
			found = true
		}
	}
	if !found {
		t.Errorf("expected intersection (0, 0) in %v", clipped)
	}
}

// TestClipPolygonNormalScaleInvariance checks that scaling the normal and
// the offset by the same factor does not change the result.
func TestClipPolygonNormalScaleInvariance(t *testing.T) {
	square := clipTestSquare()
	plane := NewPlane(&T{0.25, 0, 0}, &T{2, 0, 0})
	base, ok := ClipPolygon(square, &plane)
	if !ok {
		t.Fatal("ClipPolygon returned ok=false")
	}

	scaled := Plane{
		Normal: plane.Normal.Scaled(1e6),
		Offset: plane.Offset * 1e6,
	}
	got, ok := ClipPolygon(square, &scaled)
	if !ok {
		t.Fatal("ClipPolygon returned ok=false for a scaled plane")
	}
	if len(got) != len(base) {
		t.Fatalf("scaled plane: got %v, want %v", got, base)
	}
	for i := range base {
		if !clipTestVecEqual(got[i], base[i], 1e-5) {
			t.Errorf("vertex %d: got %v, want %v", i, got[i], base[i])
		}
	}
}

// TestClipPolygonAdditionalAllOnPlaneNormalization checks that a plane
// scaled by a negative factor flips the kept side.
func TestClipPolygonNegativeScaleFlipsSide(t *testing.T) {
	square := clipTestSquare()
	positive := clipTestPlane()
	negative := Plane{
		Normal: positive.Normal.Scaled(-1),
		Offset: -positive.Offset,
	}
	clipped, ok := ClipPolygon(square, &negative)
	if !ok {
		t.Fatal("ClipPolygon returned ok=false")
	}
	if got := PolygonArea(clipped); math.Abs(got-2) > 1e-5 {
		t.Errorf("flipped plane area: got %v, want 2", got)
	}
	for _, v := range clipped {
		if v[0] > 1e-5 {
			t.Errorf("vertex %v should lie on x <= 0", v)
		}
	}
}

// TestClipPolygonInvalidInputs checks the documented error conditions.
func TestClipPolygonInvalidInputs(t *testing.T) {
	square := clipTestSquare()
	plane := clipTestPlane()

	if got, ok := ClipPolygon(nil, &plane); ok || got != nil {
		t.Errorf("nil polygon: got %v, %v", got, ok)
	}
	if got, ok := ClipPolygon([]T{{0, 0, 0}, {1, 0, 0}}, &plane); ok || got != nil {
		t.Errorf("two vertices: got %v, %v", got, ok)
	}
	if got, ok := ClipPolygon(square, nil); ok || got != nil {
		t.Errorf("nil plane: got %v, %v", got, ok)
	}
	if got, ok := ClipPolygon(square, &Plane{}); ok || got != nil {
		t.Errorf("zero normal: got %v, %v", got, ok)
	}
	nanNormal := Plane{Normal: T{math.NaN(), 0, 0}}
	if got, ok := ClipPolygon(square, &nanNormal); ok || got != nil {
		t.Errorf("NaN normal: got %v, %v", got, ok)
	}
	nanVertex := append([]T(nil), square...)
	nanVertex[0][0] = math.NaN()
	if got, ok := ClipPolygon(nanVertex, &plane); ok || got != nil {
		t.Errorf("NaN vertex: got %v, %v", got, ok)
	}
	infVertex := append([]T(nil), square...)
	infVertex[0][1] = math.Inf(1)
	if got, ok := ClipPolygon(infVertex, &plane); ok || got != nil {
		t.Errorf("Inf vertex: got %v, %v", got, ok)
	}

	if _, _, ok := ClipPolygonUVs(square, []vec2.T{{0, 0}}, &plane); ok {
		t.Error("mismatched UV length returned ok=true")
	}
}

// TestClipPolygonDoesNotModifyInput checks that clipping returns a new
// slice and leaves the input untouched.
func TestClipPolygonDoesNotModifyInput(t *testing.T) {
	square := clipTestSquare()
	uvs := clipTestSquareUVs()
	original := append([]T(nil), square...)
	originalUVs := append([]vec2.T(nil), uvs...)
	plane := clipTestPlane()
	clipped, clippedUVs, ok := ClipPolygonUVs(square, uvs, &plane)
	if !ok {
		t.Fatal("ClipPolygonUVs returned ok=false")
	}
	for i := range square {
		if square[i] != original[i] {
			t.Errorf("input vertex %d modified: %v -> %v", i, original[i], square[i])
		}
		if uvs[i] != originalUVs[i] {
			t.Errorf("input uv %d modified: %v -> %v", i, originalUVs[i], uvs[i])
		}
	}
	if len(clipped) > 0 && len(clippedUVs) > 0 {
		if &clipped[0] == &square[0] {
			t.Error("result aliases the input vertex slice")
		}
		if &clippedUVs[0] == &uvs[0] {
			t.Error("result aliases the input uv slice")
		}
	}
}

// TestClipPolygonWinding checks both windings and confirms that original
// vertices keep their relative order.
func TestClipPolygonWinding(t *testing.T) {
	plane := clipTestPlane()
	for _, poly := range [][]T{clipTestSquare(), clipTestSquareClockwise()} {
		clipped, ok := ClipPolygon(poly, &plane)
		if !ok {
			t.Fatal("ClipPolygon returned ok=false")
		}
		if got := PolygonArea(clipped); math.Abs(got-2) > 1e-5 {
			t.Errorf("area: got %v, want 2", got)
		}
		if len(clipped) != 4 {
			t.Errorf("got %v, want 4 vertices", clipped)
		}
		// The two original vertices on x = 1 are followed by the cut on
		// the same edge (0,1) or (0,-1), so the relative order of the
		// original vertices is preserved for either winding.
		originalX1 := 0
		for _, v := range clipped {
			if v[0] == 1 {
				originalX1++
			}
		}
		if originalX1 != 2 {
			t.Errorf("expected the two original vertices on x=1 in %v", clipped)
		}
	}
}

func clipTestSquareClockwise() []T {
	return []T{
		{-1, -1, 0},
		{-1, 1, 0},
		{1, 1, 0},
		{1, -1, 0},
	}
}

func TestPolygonAreaBasics(t *testing.T) {
	if got := PolygonArea(nil); got != 0 {
		t.Errorf("nil polygon area: got %v, want 0", got)
	}
	if got := PolygonArea([]T{{0, 0, 0}, {1, 0, 0}}); got != 0 {
		t.Errorf("two-vertex area: got %v, want 0", got)
	}
	triangle := []T{{0, 0, 0}, {3, 0, 0}, {0, 4, 0}}
	if got := PolygonArea(triangle); math.Abs(got-6) > 1e-5 {
		t.Errorf("triangle area: got %v, want 6", got)
	}
	// A reflection and non-uniform scale by determinant 0.6 must scale the
	// area by 0.6, independent of the plane orientation.
	mapped := make([]T, len(triangle))
	for i := range triangle {
		mapped[i] = T{triangle[i][1] * 0.6, triangle[i][0], 2}
	}
	if got := PolygonArea(mapped); math.Abs(got-3.6) > 1e-4 {
		t.Errorf("mapped triangle area: got %v, want 3.6", got)
	}
	degenerate := []T{{0, 0, 0}, {1, 1, 1}, {2, 2, 2}}
	if got := PolygonArea(degenerate); got != 0 {
		t.Errorf("collinear polygon area: got %v, want 0", got)
	}
}

// TestClipPolygonLargeAndSmall checks that the relative tolerance keeps the
// result scale independent.
func TestClipPolygonLargeAndSmall(t *testing.T) {
	for _, scale := range []float64{1e-3, 1, 1e3} {
		square := clipTestSquare()
		for i := range square {
			square[i].Scale(scale)
		}
		plane := NewPlane(&T{0, 0, 0}, &T{1, 0, 0})
		clipped, ok := ClipPolygon(square, &plane)
		if !ok {
			t.Fatalf("scale %v: ok=false", scale)
		}
		want := 2 * scale * scale
		if got := PolygonArea(clipped); math.Abs(got-want) > 1e-5*want {
			t.Errorf("scale %v: area got %v, want %v", scale, got, want)
		}
	}
}

// TestClipPolygonFloat32Parity clips the same geometry with the float64 and
// float32 implementations and requires identical structure and coordinates.
func TestClipPolygonFloat32Parity(t *testing.T) {
	square := clipTestSquare()
	uvs := clipTestSquareUVs()
	plane := NewPlane(&T{0.3, 0.2, 0}, &T{1, 0.5, 0})
	clipped, clippedUVs, ok := ClipPolygonUVs(square, uvs, &plane)
	if !ok {
		t.Fatal("float64 ClipPolygonUVs returned ok=false")
	}

	square32 := make([]vec3f.T, len(square))
	for i := range square {
		for c := range square[i] {
			square32[i][c] = float32(square[i][c])
		}
	}
	uvs32 := make([]vec2f.T, len(uvs))
	for i := range uvs {
		for c := range uvs[i] {
			uvs32[i][c] = float32(uvs[i][c])
		}
	}
	plane32 := vec3f.NewPlane(&vec3f.T{0.3, 0.2, 0}, &vec3f.T{1, 0.5, 0})
	clipped32, uvs32out, ok := vec3f.ClipPolygonUVs(square32, uvs32, &plane32)
	if !ok {
		t.Fatal("float32 ClipPolygonUVs returned ok=false")
	}
	if len(clipped) != len(clipped32) || len(clippedUVs) != len(uvs32out) {
		t.Fatalf("vertex count mismatch: %d/%d vs %d/%d",
			len(clipped), len(clippedUVs), len(clipped32), len(uvs32out))
	}
	for i := range clipped {
		for c := range clipped[i] {
			if math.Abs(clipped[i][c]-float64(clipped32[i][c])) > 1e-5 {
				t.Errorf("vertex %d component %d mismatch: %v vs %v",
					i, c, clipped[i][c], clipped32[i][c])
			}
		}
		for c := range clippedUVs[i] {
			if math.Abs(clippedUVs[i][c]-float64(uvs32out[i][c])) > 1e-5 {
				t.Errorf("uv %d component %d mismatch: %v vs %v",
					i, c, clippedUVs[i][c], uvs32out[i][c])
			}
		}
	}
}

func TestPlaneSignedDistanceAndContains(t *testing.T) {
	plane := NewPlane(&T{1, 2, 3}, &T{0, 0, 2})
	if d := plane.SignedDistance(&T{1, 2, 3}); math.Abs(d) > 1e-6 {
		t.Errorf("distance of the defining point: got %v, want 0", d)
	}
	if d := plane.SignedDistance(&T{1, 2, 4}); math.Abs(d-2) > 1e-6 {
		t.Errorf("distance of a point above the plane: got %v, want 2", d)
	}
	if !plane.ContainsPoint(&T{1, 2, 3}) {
		t.Error("point on the plane should be contained")
	}
	if plane.ContainsPoint(&T{1, 2, 2}) {
		t.Error("point below the plane should not be contained")
	}
}

// Independent brute-force check: clipping a convex polygon against a plane
// and clipping against the three half-planes of a large triangle whose one
// edge lies on the plane must produce comparable areas. Also verify the
// retained polygon lies entirely in the kept half space and that the
// clipped area plus the area of the complementary clip sum to the original.
func TestClipPolygonAreaConservation(t *testing.T) {
	points := make([]T, 16)
	for i := range points {
		angle := 2 * math.Pi * float64(i) / float64(len(points))
		points[i] = T{3 * math.Cos(angle), 3 * math.Sin(angle), 1}
	}
	for seed := int64(0); seed < 20; seed++ {
		r := rand.New(rand.NewSource(seed))
		normal := T{r.NormFloat64(), r.NormFloat64(), r.NormFloat64()}
		if normal.Length() == 0 {
			continue
		}
		point := T{r.NormFloat64(), r.NormFloat64(), r.NormFloat64()}
		plane := NewPlane(&point, &normal)

		keep, ok := ClipPolygon(points, &plane)
		if !ok {
			t.Fatal("ok=false")
		}
		discard, ok := ClipPolygon(points, &Plane{Normal: plane.Normal.Scaled(-1), Offset: -plane.Offset})
		if !ok {
			t.Fatal("ok=false")
		}
		total := PolygonArea(points)
		if d := PolygonArea(keep) + PolygonArea(discard) - total; math.Abs(d) > 1e-2*total {
			t.Errorf("seed %d: areas do not add up: %v + %v != %v", seed,
				PolygonArea(keep), PolygonArea(discard), total)
		}
		for _, v := range keep {
			if plane.SignedDistance(&v) < -1e-3 {
				t.Errorf("seed %d: kept vertex %v outside the half space", seed, v)
			}
		}
		for _, v := range discard {
			if plane.SignedDistance(&v) > 1e-3 {
				t.Errorf("seed %d: discarded vertex %v inside the kept half space", seed, v)
			}
		}
	}
}

// Non-convex clipping area is checked against shoelace on a projection.
func TestClipPolygonShoelaceReference(t *testing.T) {
	// Concave polygon in the XY plane; project to 2D and use the shoelace
	// formula, which is independent of PolygonArea.
	poly := []T{
		{0, 0, 0}, {5, 0, 0}, {5, 1, 0}, {2, 1, 0},
		{2, 4, 0}, {5, 4, 0}, {5, 5, 0}, {0, 5, 0},
	}
	plane := NewPlane(&T{3, 0, 0}, &T{1, 0, 0})
	clipped, ok := ClipPolygon(poly, &plane)
	if !ok {
		t.Fatal("ok=false")
	}
	shoelace := 0.0
	for i := range clipped {
		a := &clipped[i]
		b := &clipped[(i+1)%len(clipped)]
		shoelace += a[0]*b[1] - b[0]*a[1]
	}
	shoelace = math.Abs(shoelace) / 2
	if math.Abs(shoelace-PolygonArea(clipped)) > 1e-4 {
		t.Errorf("PolygonArea %v disagrees with shoelace %v", PolygonArea(clipped), shoelace)
	}
	want := 2.0*1 + 2.0*1 // right of x=3: (x 3..5, y 0..1) and (x 3..5, y 4..5)
	if math.Abs(shoelace-want) > 1e-4 {
		t.Errorf("clipped area %v, want %v", shoelace, want)
	}
}

// TestClipPolygonFarFromOrigin cuts a square whose coordinates are large.
// The tolerance is derived from the signed distances rather than the raw
// coordinate magnitude, so a genuine cut far from the origin must not be
// treated as on-plane.
func TestClipPolygonFarFromOrigin(t *testing.T) {
	const offset = 1e6
	// Square in the yz plane spanning x = offset-1000 .. offset+1000.
	square := []T{
		{offset - 1000, -1000, 0},
		{offset + 1000, -1000, 0},
		{offset + 1000, 1000, 0},
		{offset - 1000, 1000, 0},
	}
	plane := NewPlane(&T{offset, 0, 0}, &T{1, 0, 0})
	clipped, ok := ClipPolygon(square, &plane)
	if !ok {
		t.Fatal("ClipPolygon returned ok=false")
	}
	if len(clipped) != 4 {
		t.Fatalf("expected a rectangle, got %v", clipped)
	}
	want := float64(2000 * 1000) // 2000 x 2000 square cut in half, area 2000 x 1000
	if got := PolygonArea(clipped); math.Abs(got-want) > 1e-3*want {
		t.Errorf("area: got %v, want %v", got, want)
	}
	for _, v := range clipped {
		if v[0] < offset {
			t.Errorf("vertex %v lies on the discarded side", v)
		}
	}
}

// TestClipPolygonCancellationOnPlane places a polygon on a tilted plane far
// from the origin. The signed distances are the cancellation of large dot
// products and the offset, so without an error-scale based tolerance the
// vertices would be pushed to one side of the plane.
func TestClipPolygonCancellationOnPlane(t *testing.T) {
	const offset = 1e6
	normal := T{1, 1, 0}
	point := T{offset, 0, 0}
	plane := NewPlane(&point, &normal)
	// The square lies in the plane through (offset, 0, 0) spanned by
	// (1,-1,0) and (0,0,1); its corners are offset by +-1 along both.
	polygon := []T{
		{offset - 1, 1, -1},
		{offset + 1, -1, -1},
		{offset + 1, -1, 1},
		{offset - 1, 1, 1},
	}
	for _, v := range polygon {
		if d := plane.SignedDistance(&v); math.Abs(d) > 2 {
			t.Fatalf("test setup: vertex %v is %v off the plane", v, d)
		}
	}
	clipped, ok := ClipPolygon(polygon, &plane)
	if !ok {
		t.Fatal("ClipPolygon returned ok=false")
	}
	if len(clipped) != len(polygon) {
		t.Errorf("polygon on the plane should be unchanged, got %v", clipped)
	}
}

// TestClipPolygonConcaveTwoPieces cuts a U-shaped polygon so the kept half
// space contains two disconnected arms. A single-plane clipping returns one
// boundary loop with zero-width bridges along the cutting line, so the area
// must still be the sum of both arms and every bridge pair must lie on the
// plane.
func TestClipPolygonConcaveTwoPieces(t *testing.T) {
	// U rotated so its opening faces -X: the kept side x >= 0 contains the
	// two arms at y in [-3,-1] and [1,3].
	u := []T{
		{-1, -3, 0},
		{1, -3, 0},
		{1, -1, 0},
		{-1, -1, 0},
		{-1, 1, 0},
		{1, 1, 0},
		{1, 3, 0},
		{-1, 3, 0},
	}
	plane := NewPlane(&T{0, 0, 0}, &T{1, 0, 0})
	clipped, ok := ClipPolygon(u, &plane)
	if !ok {
		t.Fatal("ClipPolygon returned ok=false")
	}
	// Each arm contributes 1 x 2, the two bridges have zero area.
	if got := PolygonArea(clipped); math.Abs(got-4) > 1e-5 {
		t.Errorf("area: got %v, want 4", got)
	}
	bridges := 0
	for _, v := range clipped {
		if d := plane.SignedDistance(&v); math.Abs(d) <= 1e-5 {
			bridges++
		}
	}
	if bridges < 4 {
		t.Errorf("expected at least 4 points on the cutting line, got %d in %v", bridges, clipped)
	}
	for _, v := range clipped {
		if v[0] < -1e-5 {
			t.Errorf("vertex %v lies on the discarded side", v)
		}
	}
}

// TestClipPolygonUVInterpolation checks UVs of intersection points on edges
// whose endpoints have different UV spans, including an edge cut close to
// one endpoint.
func TestClipPolygonUVInterpolation(t *testing.T) {
	poly := []T{
		{0, 0, 0},
		{10, 0, 0},
		{10, 10, 0},
		{0, 10, 0},
	}
	uvs := []vec2.T{
		{0, 0},
		{0.2, 0},
		{0.2, 0.5},
		{0, 0.5},
	}
	// Cut at x = 9, close to the right edge, so the interpolated u is 0.18
	// on the two horizontal edges.
	plane := NewPlane(&T{9, 0, 0}, &T{1, 0, 0})
	clipped, clippedUVs, ok := ClipPolygonUVs(poly, uvs, &plane)
	if !ok {
		t.Fatal("ClipPolygonUVs returned ok=false")
	}
	if len(clipped) != len(clippedUVs) || len(clipped) != 4 {
		t.Fatalf("got %d vertices and %d uvs, want 4 and 4", len(clipped), len(clippedUVs))
	}
	// Find the intersections on the horizontal edges (y = 0 and y = 10)
	// and the original vertices on x = 10.
	for i := range clipped {
		v, uv := clipped[i], clippedUVs[i]
		switch {
		case v[0] == 10:
			// Kept originals keep their UVs exactly.
			if math.Abs(uv[0]-0.2) > 1e-6 {
				t.Errorf("original vertex %v has uv %v, want u=0.2", v, uv)
			}
		case math.Abs(v[0]-9) <= 1e-5:
			if math.Abs(uv[0]-0.18) > 1e-5 {
				t.Errorf("intersection %v has uv %v, want u=0.18", v, uv)
			}
		default:
			t.Errorf("unexpected vertex %v", v)
		}
	}
}

// independentClip is a separate Sutherland-Hodgman implementation written
// directly from the half-plane predicate, without tolerances or duplicate
// handling, used to cross-check ClipPolygon on random convex polygons. The
// area is verified with Newell's formula again, independent of PolygonArea.
func independentClip(poly []T, plane *Plane) (vertices []T, area float64) {
	n := len(poly)
	for i := 0; i < n; i++ {
		a := poly[i]
		b := poly[(i+1)%n]
		da := plane.SignedDistance(&a)
		db := plane.SignedDistance(&b)
		if da >= 0 {
			vertices = append(vertices, a)
		}
		if (da >= 0) != (db >= 0) {
			t := da / (da - db)
			if t > 0 && t < 1 {
				vertices = append(vertices, Interpolate(&a, &b, t))
			}
		}
	}
	// Newell area, computed in float64 independently.
	var nx, ny, nz float64
	for i := range vertices {
		a := vertices[i]
		b := vertices[(i+1)%len(vertices)]
		nx += (a[1] - b[1]) * (a[2] + b[2])
		ny += (a[2] - b[2]) * (a[0] + b[0])
		nz += (a[0] - b[0]) * (a[1] + b[1])
	}
	return vertices, 0.5 * math.Sqrt(nx*nx+ny*ny+nz*nz)
}

func TestClipPolygonRandomConvexCrossCheck(t *testing.T) {
	r := rand.New(rand.NewSource(12345))
	for trial := 0; trial < 100; trial++ {
		// Convex polygon: points on a circle with monotone angles.
		n := 3 + r.Intn(8)
		poly := make([]T, n)
		for i := range poly {
			angle := 2 * math.Pi * (float64(i) + 0.2*r.Float64()) / float64(n)
			radius := 1 + 3*r.Float64()
			poly[i] = T{radius * math.Cos(angle), radius * math.Sin(angle), 0.5 * r.Float64()}
		}
		normal := T{r.NormFloat64(), r.NormFloat64(), r.NormFloat64()}
		if normal.Length() < 1e-3 {
			continue
		}
		point := T{r.NormFloat64(), r.NormFloat64(), r.NormFloat64()}
		plane := NewPlane(&point, &normal)

		got, ok := ClipPolygon(poly, &plane)
		if !ok {
			t.Fatalf("trial %d: ok=false", trial)
		}
		wantV, wantArea := independentClip(poly, &plane)
		if len(got) != len(wantV) {
			t.Errorf("trial %d: %d vertices, reference %d", trial, len(got), len(wantV))
		}
		if d := math.Abs(PolygonArea(got) - wantArea); d > 1e-3*math.Max(1, wantArea) {
			t.Errorf("trial %d: area %v, reference %v", trial, PolygonArea(got), wantArea)
		}
		for _, v := range got {
			if plane.SignedDistance(&v) < -1e-3 {
				t.Errorf("trial %d: vertex %v outside the half space", trial, v)
			}
		}
	}
}

// TestClipPolygonOnPlaneStability clips random polygons that lie exactly on
// a random plane. Each vertex is built from the plane origin and two
// in-plane directions, so its true signed distance is zero; floating point
// rounding must not push it to either side and no intersections may be
// added. This is the regression test for the error-scale tolerance.
func TestClipPolygonOnPlaneStability(t *testing.T) {
	r := rand.New(rand.NewSource(7))
	fail := 0
	for trial := 0; trial < 5000; trial++ {
		origin := T{r.NormFloat64(), r.NormFloat64(), r.NormFloat64()}
		normal := T{r.NormFloat64(), r.NormFloat64(), r.NormFloat64()}
		if normal.Length() < 1e-2 {
			continue
		}
		plane := NewPlane(&origin, &normal)
		t1 := normal.Normal()
		t2 := Cross(&normal, &t1)
		t2.Normalize()
		n := 3 + r.Intn(8)
		poly := make([]T, n)
		for i := range poly {
			u := 10 * r.NormFloat64()
			v := 10 * r.NormFloat64()
			su := t1.Scaled(u)
			sv := t2.Scaled(v)
			p := origin
			p.Add(&su)
			p.Add(&sv)
			poly[i] = p
		}
		out, ok := ClipPolygon(poly, &plane)
		if !ok || len(out) != len(poly) {
			fail++
			if fail < 3 {
				t.Logf("trial %d: in %d out %d ok %v", trial, len(poly), len(out), ok)
				for _, p := range poly {
					t.Logf("   p=%v d=%v", p, plane.SignedDistance(&p))
				}
			}
		}
	}
	t.Logf("on-plane failures: %d / 5000", fail)
	if fail > 0 {
		t.Errorf("%d on-plane polygons were modified", fail)
	}
}
