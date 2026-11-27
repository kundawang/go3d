package vec3

import (
	"github.com/ungerik/go3d/float64/vec2"
	"math"
)

// clipPolygonEpsilon is the relative tolerance used to decide whether a
// vertex lies on the clipping plane and whether two points coincide.
//
// It is a few floating point rounding errors large. Tolerances are derived
// from the magnitude of the terms that produce a signed distance, so
// coordinates that should lie exactly on the plane are not pushed to either
// side by rounding, while a genuine offset larger than a few units in the
// last place is still classified correctly.
const clipPolygonEpsilon float64 = 1e-14

// ClipPolygon clips a polygon against a single oriented plane and returns
// the part of the polygon on the positive side of the plane.
//
// The polygon is given by its vertices in boundary order; clockwise and
// counterclockwise winding are both accepted and the winding is preserved.
// The positive (kept) side is the set of points p with
// plane.SignedDistance(p) >= 0, so vertices and intersection points exactly
// on the plane are kept.
//
// The result is a new polygon in the same boundary order as the input: kept
// original vertices keep their relative order and a new intersection point
// is inserted between the endpoints of the edge it lies on. Original
// vertices are never reordered or removed, so a polygon that lies entirely
// on the positive side is returned as a copy with the same vertex order.
//
// For a concave polygon the intersection with one half space can consist of
// several disconnected pieces. As in any single-plane Sutherland-Hodgman
// implementation, such pieces are returned as one boundary loop that is
// connected along the cutting line by zero-width bridges; the pieces can be
// told apart by the pairs of intersection points that lie on the plane.
//
// Vertices whose signed distance is within a small relative tolerance of the
// rounding error of that distance are treated as lying on the plane. The
// tolerance is derived per vertex from the terms that cancel in the distance
// computation, so it works both near and far from the origin. This keeps the
// result stable for vertices that are on the plane up to rounding, for edges
// that are almost parallel to the plane, and for polygons that only touch the
// plane at a vertex or along an edge. An intersection point that would
// coincide with an endpoint or with the previously emitted vertex is
// represented by that vertex instead of a duplicate.
//
// The boolean result is false if the input cannot be clipped, in which case
// the returned slice is nil: fewer than three vertices, a nil plane, a plane
// with a zero, NaN or infinite normal, or a vertex, normal, offset or
// signed distance that is NaN or infinite. Otherwise the result is valid
// even if it is empty (the whole polygon lies on the negative side) or has
// fewer than three vertices (the retained part degenerated to an edge or a
// point).
//
// The input slice is not modified and the result never aliases it.
func ClipPolygon(polygon []T, plane *Plane) ([]T, bool) {
	clipped, _, ok := clipPolygon(polygon, nil, plane)
	return clipped, ok
}

// ClipPolygonUVs is like ClipPolygon but also clips one per-vertex UV
// coordinate (vec2.T) per polygon vertex. The UV coordinate of a new
// intersection point is interpolated with the same edge parameter t as its
// position, so an edge is cut at one consistent point in position and
// texture space.
//
// len(uvs) must equal len(polygon); otherwise the boolean result is false.
// The returned UV slice always has the same length as the returned polygon.
// The input slices are not modified.
func ClipPolygonUVs(polygon []T, uvs []vec2.T, plane *Plane) ([]T, []vec2.T, bool) {
	if len(uvs) != len(polygon) {
		return nil, nil, false
	}
	return clipPolygon(polygon, uvs, plane)
}

// PolygonArea returns the area of a polygon whose vertices are given in
// boundary order.
//
// It works for planar polygons of arbitrary orientation in 3D and for both
// windings; the result is always non-negative. The polygon need not be
// convex. The area is computed with Newell's method, which needs no
// projection direction and stays stable for almost degenerate polygons.
// Fewer than three vertices return zero; vertices with non-finite
// coordinates return NaN.
func PolygonArea(polygon []T) float64 {
	if len(polygon) < 3 {
		return 0
	}
	var nx, ny, nz float64
	for i := range polygon {
		a := &polygon[i]
		b := &polygon[(i+1)%len(polygon)]
		nx += (a[1] - b[1]) * (a[2] + b[2])
		ny += (a[2] - b[2]) * (a[0] + b[0])
		nz += (a[0] - b[0]) * (a[1] + b[1])
	}
	return 0.5 * math.Sqrt(nx*nx+ny*ny+nz*nz)
}

// clipPolygon implements Sutherland-Hodgman clipping against a single
// plane. uvs is either nil or has the same length as polygon.
func clipPolygon(polygon []T, uvs []vec2.T, plane *Plane) ([]T, []vec2.T, bool) {
	if len(polygon) < 3 {
		return nil, nil, false
	}
	if uvs != nil && len(uvs) != len(polygon) {
		return nil, nil, false
	}
	if plane == nil {
		return nil, nil, false
	}

	// Normalize the plane so the signed distances, and with them the
	// tolerance, do not depend on the magnitude of Normal and Offset. A
	// length that is zero or not finite (NaN included) makes the plane
	// unusable.
	normalLength := plane.Normal.Length()
	if !(normalLength > 0) || !clipFiniteFloat(normalLength) {
		return nil, nil, false
	}
	normal := plane.Normal.Scaled(1 / normalLength)
	offset := plane.Offset / normalLength
	if !clipFiniteVec(&normal) || !clipFiniteFloat(offset) {
		return nil, nil, false
	}

	// Precompute the signed distances and, per vertex, the magnitude of the
	// products that produce them. A vertex on the plane is computed as the
	// cancellation of the products normal[c]*vertex[c] and offset, so its
	// rounding error scales with the sum of their magnitudes, not with the
	// resulting distance. Deriving the tolerance from that scale keeps
	// vertices on the plane stable even far from the origin, while vertices
	// genuinely off the plane still classify correctly.
	distances := make([]float64, len(polygon))
	tolerances := make([]float64, len(polygon))
	for i := range polygon {
		vertex := &polygon[i]
		if !clipFiniteVec(vertex) {
			return nil, nil, false
		}
		d := Dot(&normal, vertex)
		if !clipFiniteFloat(d) {
			return nil, nil, false
		}
		computationScale := math.Abs(offset)
		for c := range normal {
			computationScale += math.Abs(normal[c] * vertex[c])
		}
		d += offset
		if !clipFiniteFloat(d) {
			return nil, nil, false
		}
		distances[i] = d
		tolerances[i] = clipPolygonEpsilon * computationScale
	}

	hasUVs := uvs != nil
	out := make([]T, 0, len(polygon)+1)
	isIntersection := make([]bool, 0, len(polygon)+1)
	var outUVs []vec2.T
	if hasUVs {
		outUVs = make([]vec2.T, 0, len(polygon)+1)
	}

	// emitOriginal appends an original vertex. Original vertices are
	// never removed or reordered, so a polygon that lies completely on
	// the positive side is returned as a copy.
	emitOriginal := func(position T, uv vec2.T) {
		out = append(out, position)
		isIntersection = append(isIntersection, false)
		if hasUVs {
			outUVs = append(outUVs, uv)
		}
	}

	n := len(polygon)
	for i := 0; i < n; i++ {
		current := &polygon[i]
		next := &polygon[(i+1)%n]
		currentUV := vec2.Zero
		nextUV := vec2.Zero
		if hasUVs {
			currentUV = uvs[i]
			nextUV = uvs[(i+1)%n]
		}

		dCurrent := distances[i]
		dNext := distances[(i+1)%n]
		currentInside := dCurrent >= -tolerances[i]
		nextInside := dNext >= -tolerances[(i+1)%n]

		if currentInside {
			emitOriginal(*current, currentUV)
		}

		if currentInside != nextInside {
			// The endpoints are classified on opposite sides, so their
			// signed distances lie on opposite sides of the tolerance
			// band and the denominator cannot be zero. Mathematically t
			// lies in [0, 1]; the clamp only absorbs rounding.
			t := dCurrent / (dCurrent - dNext)
			if t <= 0 || t >= 1 {
				// The intersection coincides with an endpoint. If that
				// endpoint is inside it is emitted as an original vertex,
				// either already or in the next iteration.
				continue
			}
			position := Interpolate(current, next, t)
			if n := len(out); n > 0 && clipPointsCoincide(&out[n-1], &position) {
				// Coincides with the vertex emitted just before; keep the
				// original vertex instead of inserting a duplicate.
				continue
			}
			if nextInside && clipPointsCoincide(next, &position) {
				// Coincides with next, which is emitted as an original
				// vertex in the next iteration.
				continue
			}
			out = append(out, position)
			isIntersection = append(isIntersection, true)
			if hasUVs {
				outUVs = append(outUVs, vec2.Interpolate(&currentUV, &nextUV, t))
			}
		}
	}

	// The polygon boundary is cyclic, so the last vertex and the first
	// vertex may also be duplicates of each other, e.g. when the plane
	// touches the polygon at a single corner. Only an intersection point
	// is dropped; two coincident original vertices are left untouched.
	for len(out) > 1 && clipPointsCoincide(&out[0], &out[len(out)-1]) {
		switch {
		case isIntersection[len(out)-1]:
			last := len(out) - 1
			out = out[:last]
			isIntersection = isIntersection[:last]
			if hasUVs {
				outUVs = outUVs[:last]
			}
		case isIntersection[0]:
			out = out[1:]
			isIntersection = isIntersection[1:]
			if hasUVs {
				outUVs = outUVs[1:]
			}
		default:
			return out, outUVs, true
		}
	}

	return out, outUVs, true
}

// clipPointsCoincide reports whether a and b agree within a small relative
// tolerance in every component. Comparing each component relative to its
// own magnitude keeps the test meaningful for coordinates far from the
// origin that do not influence the clipping plane.
func clipPointsCoincide(a, b *T) bool {
	for c := range a {
		if math.Abs(a[c]-b[c]) > clipPolygonEpsilon*clipMagnitude(&a[c], &b[c]) {
			return false
		}
	}
	return true
}

// clipMagnitude returns the larger absolute value of a and b, the scale
// against which their difference is compared.
func clipMagnitude(a, b *float64) float64 {
	absA, absB := math.Abs(*a), math.Abs(*b)
	if absA > absB {
		return absA
	}
	return absB
}

// clipFiniteFloat reports whether f is neither NaN nor an infinity.
func clipFiniteFloat(f float64) bool {
	// NaN is the only value that does not compare equal to itself, and
	// any magnitude above math.MaxFloat64 is an infinity.
	return f == f && -math.MaxFloat64 <= f && f <= math.MaxFloat64
}

// clipFiniteVec reports whether all components of v are finite.
func clipFiniteVec(v *T) bool {
	for c := range v {
		if !clipFiniteFloat(v[c]) {
			return false
		}
	}
	return true
}
