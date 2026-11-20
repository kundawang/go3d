// Package frustum provides view frustum culling for float64 types.
//
// A frustum is extracted from a combined view-projection matrix and can
// classify axis-aligned bounding boxes (vec3.Box), optionally transformed
// by an arbitrary mat4, as Inside, Intersect or Outside the frustum.
// Only Outside boxes may be culled by a renderer.
//
// The extraction follows Gribb & Hartmann ("Fast Extraction of Viewing
// Frustum Planes from the World-View-Projection Matrix") and assumes the
// OpenGL clip space convention (-w <= x, y, z <= w) used by the mat4
// package. It works unchanged for left- and right-handed view matrices
// and for perspective as well as orthographic projections, because
// handedness is already absorbed into the combined matrix.
package frustum

import (
	"errors"
	"fmt"

	"math"
	"github.com/ungerik/go3d/float64/mat4"
	"github.com/ungerik/go3d/float64/vec3"
	"github.com/ungerik/go3d/float64/vec4"
)

// Relation describes the spatial relation of an object to a frustum.
type Relation int

const (
	// Outside means the object lies completely outside the frustum.
	// Only objects classified Outside may be culled.
	Outside Relation = iota
	// Intersect means the object intersects at least one frustum plane,
	// so it is partially (or, conservatively, possibly) visible.
	// Objects exactly tangent to a plane are classified Intersect,
	// never Outside.
	Intersect
	// Inside means the object lies completely inside the frustum.
	Inside
)

// String formats Relation as string.
func (r Relation) String() string {
	switch r {
	case Outside:
		return "Outside"
	case Intersect:
		return "Intersect"
	case Inside:
		return "Inside"
	default:
		return fmt.Sprintf("Relation(%d)", int(r))
	}
}

// Indices of the frustum planes in T.Planes.
const (
	PlaneLeft = iota
	PlaneRight
	PlaneBottom
	PlaneTop
	PlaneNear
	PlaneFar
	// PlaneCount is the number of frustum planes.
	PlaneCount
)

var (
	// ErrNonFiniteMatrix is returned by FromViewProjection when the
	// view-projection matrix contains NaN or infinite elements. This
	// covers degenerate projection parameters such as znear == zfar,
	// which produce infinities in the projection matrix.
	ErrNonFiniteMatrix = errors.New("frustum: view-projection matrix contains non-finite values")
	// ErrDegenerateMatrix is returned by FromViewProjection when the
	// matrix yields a plane with a (near) zero normal, e.g. a zero
	// matrix or one whose clip-space volume collapsed.
	ErrDegenerateMatrix = errors.New("frustum: view-projection matrix produces a degenerate plane")
)

// T represents a view frustum by its 6 bounding planes.
//
// Each plane is stored as (a, b, c, d) with a*x+b*y+c*z+d >= 0 on the
// inside of the frustum. Normals (a, b, c) point inwards and are
// normalized to unit length, so the plane equation evaluates to the
// signed distance from the plane.
type T struct {
	Planes [PlaneCount]vec4.T
}

// FromViewProjection extracts the frustum of a combined
// view-projection matrix (projection * view).
//
// The matrix may be left- or right-handed and may contain a perspective
// or orthographic projection, as long as it maps to OpenGL-style clip
// space (-w <= x, y, z <= w), which is what the mat4 package's
// AssignPerspective, AssignFrustum and AssignOrthogonalProjection
// produce.
//
// An error is returned for matrices containing NaN or infinite elements
// (ErrNonFiniteMatrix, e.g. from znear == zfar) and for matrices that
// produce a degenerate plane (ErrDegenerateMatrix, e.g. the zero
// matrix). See also MustFromViewProjection.
func FromViewProjection(vp *mat4.T) (T, error) {
	var f T
	for c := 0; c < 4; c++ {
		for r := 0; r < 4; r++ {
			if v := vp[c][r]; math.IsNaN(v) || math.IsInf(v, 0) {
				return f, ErrNonFiniteMatrix
			}
		}
	}

	// Rows of the column-major matrix.
	row0 := vec4.T{vp[0][0], vp[1][0], vp[2][0], vp[3][0]}
	row1 := vec4.T{vp[0][1], vp[1][1], vp[2][1], vp[3][1]}
	row2 := vec4.T{vp[0][2], vp[1][2], vp[2][2], vp[3][2]}
	row3 := vec4.T{vp[0][3], vp[1][3], vp[2][3], vp[3][3]}

	planes := [PlaneCount]vec4.T{
		PlaneLeft:   add4(row3, row0),
		PlaneRight:  sub4(row3, row0),
		PlaneBottom: add4(row3, row1),
		PlaneTop:    sub4(row3, row1),
		PlaneNear:   add4(row3, row2),
		PlaneFar:    sub4(row3, row2),
	}
	for i := range planes {
		p := &planes[i]
		length := math.Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
		// !(length > eps) also catches NaN; Inf would normalize to NaN.
		if !(length > 1e-280) || math.IsInf(length, 0) {
			return f, fmt.Errorf("%w (plane %d)", ErrDegenerateMatrix, i)
		}
		inv := 1 / length
		f.Planes[i] = vec4.T{p[0] * inv, p[1] * inv, p[2] * inv, p[3] * inv}
	}
	return f, nil
}

// MustFromViewProjection is like FromViewProjection but panics on error.
// Useful for frustums built from known-valid matrices, e.g. in tests or
// during initialization.
func MustFromViewProjection(vp *mat4.T) T {
	f, err := FromViewProjection(vp)
	if err != nil {
		panic(err)
	}
	return f
}

// ClassifyBox classifies an axis-aligned box against the frustum.
// It is equivalent to ClassifyTransformedBox with the identity matrix.
func (f *T) ClassifyBox(box *vec3.Box) Relation {
	lo, hi, ok := orderedFiniteBox(box)
	if !ok {
		// Conservative: geometry that cannot be classified is never culled.
		return Intersect
	}
	relation := Inside
	for i := 0; i < PlaneCount; i++ {
		relation = f.classifyPlane(&f.Planes[i], &lo, &hi, relation)
		if relation == Outside {
			return Outside
		}
	}
	return relation
}

// ClassifyTransformedBox classifies an axis-aligned box that is
// transformed by the given matrix (e.g. a model matrix) against the
// frustum, without computing the transformed corners.
//
// The transform may contain rotation, translation, uniform or
// non-uniform scaling (including negative scaling, i.e. mirroring) and
// shear. Singular affine transforms (e.g. a zero scale axis collapsing
// the box to a plane, line or point) are handled exactly. The transform
// is expected to be affine, i.e. its last row must be (0, 0, 0, 1),
// which holds for all typical model matrices.
//
// If the box or the transform contains NaN or infinite values, the
// result is Intersect: unclassifiable geometry is conservatively kept
// instead of being culled.
func (f *T) ClassifyTransformedBox(box *vec3.Box, transform *mat4.T) Relation {
	lo, hi, ok := orderedFiniteBox(box)
	if !ok || !finiteMat4(transform) {
		// Conservative: geometry that cannot be classified is never culled.
		return Intersect
	}
	relation := Inside
	for i := 0; i < PlaneCount; i++ {
		// Transform the world-space plane into the local space of the
		// box instead of transforming the 8 corners into world space:
		// for p_world = M * p_local the plane equation gives
		// dot(w, p_world) = dot(w * M, p_local), so the local plane is
		// the world plane (as row vector) times the matrix, i.e. its
		// component c is the dot product with matrix column c.
		w := &f.Planes[i]
		local := vec4.T{
			w[0]*transform[0][0] + w[1]*transform[0][1] + w[2]*transform[0][2] + w[3]*transform[0][3],
			w[0]*transform[1][0] + w[1]*transform[1][1] + w[2]*transform[1][2] + w[3]*transform[1][3],
			w[0]*transform[2][0] + w[1]*transform[2][1] + w[2]*transform[2][2] + w[3]*transform[2][3],
			w[0]*transform[3][0] + w[1]*transform[3][1] + w[2]*transform[3][2] + w[3]*transform[3][3],
		}
		relation = f.classifyPlane(&local, &lo, &hi, relation)
		if relation == Outside {
			return Outside
		}
	}
	return relation
}

// ContainsPoint returns whether a point lies inside or on the frustum.
// Points exactly on a plane are considered inside. A point with
// non-finite coordinates is never contained.
func (f *T) ContainsPoint(p *vec3.T) bool {
	for c := 0; c < 3; c++ {
		if math.IsNaN(p[c]) || math.IsInf(p[c], 0) {
			return false
		}
	}
	for i := 0; i < PlaneCount; i++ {
		plane := &f.Planes[i]
		if plane[0]*p[0]+plane[1]*p[1]+plane[2]*p[2]+plane[3] < 0 {
			return false
		}
	}
	return true
}

// classifyPlane updates the relation of a box (lo, hi) for a single
// plane using the positive/negative vertex test. The outside comparison
// is strict (< 0) so boxes exactly tangent to a plane are never culled.
func (f *T) classifyPlane(plane *vec4.T, lo, hi *vec3.T, relation Relation) Relation {
	// Positive vertex: the corner farthest along the plane normal.
	// Negative vertex: the opposite corner.
	var pv, nv vec3.T
	for c := 0; c < 3; c++ {
		if plane[c] >= 0 {
			pv[c] = hi[c]
			nv[c] = lo[c]
		} else {
			pv[c] = lo[c]
			nv[c] = hi[c]
		}
	}
	if plane[0]*pv[0]+plane[1]*pv[1]+plane[2]*pv[2]+plane[3] < 0 {
		return Outside
	}
	// The negative-vertex test is non-strict so that boxes exactly
	// tangent to the plane are classified Intersect, not Inside.
	if plane[0]*nv[0]+plane[1]*nv[1]+plane[2]*nv[2]+plane[3] <= 0 {
		return Intersect
	}
	return relation
}

// orderedFiniteBox returns the component-wise ordered corners of the
// box and whether all coordinates are finite. Ordering tolerates boxes
// with Min > Max components; zero-volume boxes (Min == Max, e.g. points)
// are valid and classify like any other box.
func orderedFiniteBox(box *vec3.Box) (lo, hi vec3.T, ok bool) {
	for c := 0; c < 3; c++ {
		min, max := box.Min[c], box.Max[c]
		if math.IsNaN(min) || math.IsInf(min, 0) || math.IsNaN(max) || math.IsInf(max, 0) {
			return lo, hi, false
		}
		if min <= max {
			lo[c], hi[c] = min, max
		} else {
			lo[c], hi[c] = max, min
		}
	}
	return lo, hi, true
}

func finiteMat4(m *mat4.T) bool {
	for c := 0; c < 4; c++ {
		for r := 0; r < 4; r++ {
			if v := m[c][r]; math.IsNaN(v) || math.IsInf(v, 0) {
				return false
			}
		}
	}
	return true
}

func add4(a, b vec4.T) vec4.T {
	return vec4.T{a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]}
}

func sub4(a, b vec4.T) vec4.T {
	return vec4.T{a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]}
}
