package mat4

import (
	math "github.com/ungerik/go3d/fmath"
	"github.com/ungerik/go3d/vec3"
	"github.com/ungerik/go3d/vec4"
)

const (
	// minClipW is the smallest absolute homogeneous clip-space w that Project()
	// still divides by. A clip point with |w| below this threshold lies on (or
	// practically on) the camera plane, where perspective division is undefined.
	minClipW float32 = 1e-7

	// invertSingularThreshold scales the pivot test of invert4x4():
	// a pivot smaller than this fraction of the largest absolute matrix
	// element marks the matrix as singular.
	invertSingularThreshold float32 = 1e-9

	// clipTolerance relaxes the normalized device coordinate cube test,
	// so that points landing exactly on a clipping plane are not reported
	// as being outside of the clip volume due to rounding in the projection
	// matrices themselves.
	clipTolerance float32 = 1e-5
)

// Viewport describes a non-zero-origin editor viewport in top-left pixel
// coordinates. X and Y are the pixel coordinates of the top-left corner,
// Width and Height the dimensions in pixels. Width and Height must both be
// greater than zero for Project and Unproject to succeed.
type Viewport struct {
	X, Y, Width, Height float32
}

// finiteFloat returns true if f is neither NaN nor an infinity.
func finiteFloat(f float32) bool {
	// NaN is the only value that does not compare equal to itself,
	// and any value with a magnitude above MaxFloat32 is infinite.
	return f == f && -math.MaxFloat32 <= f && f <= math.MaxFloat32
}

// finiteMat returns true if all 16 elements of mat are finite.
func finiteMat(mat *T) bool {
	for col := 0; col < 4; col++ {
		for row := 0; row < 4; row++ {
			if !finiteFloat(mat[col][row]) {
				return false
			}
		}
	}
	return true
}

// outsideClipCube reports whether a normalized device coordinate lies outside
// of the [-1, 1] cube (the clip volume after perspective division).
func outsideClipCube(ndcX, ndcY, ndcZ float32) bool {
	return ndcX < -1-clipTolerance || ndcX > 1+clipTolerance ||
		ndcY < -1-clipTolerance || ndcY > 1+clipTolerance ||
		ndcZ < -1-clipTolerance || ndcZ > 1+clipTolerance
}

// invert4x4 inverts mat using Gauss-Jordan elimination with partial pivoting.
// It returns the inverted matrix and ok == true for an invertible matrix.
// For a singular matrix (or non-finite results) it returns the zero matrix
// and ok == false.
func invert4x4(mat *T) (inv T, ok bool) {
	// Augmented matrix [mat | ident], stored row-major.
	var a [4][8]float32
	maxAbs := float32(0)
	for row := 0; row < 4; row++ {
		for col := 0; col < 4; col++ {
			v := mat[col][row]
			a[row][col] = v
			if abs := math.Abs(v); abs > maxAbs {
				maxAbs = abs
			}
		}
		a[row][row+4] = 1
	}
	if maxAbs == 0 || !finiteFloat(maxAbs) {
		return Zero, false
	}
	pivotEps := invertSingularThreshold * maxAbs

	for i := 0; i < 4; i++ {
		// Partial pivoting: find the row with the largest absolute pivot.
		pivotRow := i
		pivotAbs := math.Abs(a[i][i])
		for row := i + 1; row < 4; row++ {
			if abs := math.Abs(a[row][i]); abs > pivotAbs {
				pivotAbs = abs
				pivotRow = row
			}
		}
		if pivotAbs <= pivotEps {
			return Zero, false
		}
		if pivotRow != i {
			a[i], a[pivotRow] = a[pivotRow], a[i]
		}

		pivot := a[i][i]
		for j := 0; j < 8; j++ {
			a[i][j] /= pivot
		}
		for row := 0; row < 4; row++ {
			if row == i {
				continue
			}
			factor := a[row][i]
			for j := 0; j < 8; j++ {
				a[row][j] -= factor * a[i][j]
			}
		}
	}

	for row := 0; row < 4; row++ {
		for col := 0; col < 4; col++ {
			v := a[row][col+4]
			if !finiteFloat(v) {
				return Zero, false
			}
			inv[col][row] = v
		}
	}
	return inv, true
}

// Project transforms the world-space point world into editor viewport
// coordinates using the given view and projection matrices and viewport.
//
// The returned screen vector holds top-left origin pixel coordinates in its
// X and Y components and the window-space depth in its Z component, which is
// in the range [0, 1] for points inside the clip volume (0 on the near
// clipping plane, 1 on the far clipping plane).
//
// ok is true when the point is inside the clip volume and the projection is
// well defined. ok is false when the point is outside of the clip volume,
// when the homogeneous clip-space w is nearly zero (the point lies on the
// camera plane), when the viewport has a zero or negative width or height, or
// when any input is non-finite. In every failure case the returned vector is
// vec3.Zero.
//
// float32 and float64 implementations of this function have identical
// semantics, only their precision differs.
func Project(world *vec3.T, view, projection *T, viewport Viewport) (screen vec3.T, ok bool) {
	if !finiteFloat(world[0]) || !finiteFloat(world[1]) || !finiteFloat(world[2]) ||
		!finiteMat(view) || !finiteMat(projection) ||
		!finiteFloat(viewport.X) || !finiteFloat(viewport.Y) ||
		!finiteFloat(viewport.Width) || !finiteFloat(viewport.Height) {
		return vec3.Zero, false
	}
	if viewport.Width <= 0 || viewport.Height <= 0 {
		return vec3.Zero, false
	}

	var viewProjection T
	viewProjection.AssignMul(projection, view)

	clip := viewProjection.MulVec4(&vec4.T{world[0], world[1], world[2], 1})
	if math.Abs(clip[3]) < minClipW {
		return vec3.Zero, false
	}

	// Perspective division yields normalized device coordinates in [-1, 1].
	ndcX := clip[0] / clip[3]
	ndcY := clip[1] / clip[3]
	ndcZ := clip[2] / clip[3]
	if !finiteFloat(ndcX) || !finiteFloat(ndcY) || !finiteFloat(ndcZ) {
		return vec3.Zero, false
	}

	// Map NDC to top-left origin pixel coordinates and [0, 1] depth.
	// The Y axis is flipped because NDC +Y points up while pixel +Y points down.
	if outsideClipCube(ndcX, ndcY, ndcZ) {
		return vec3.Zero, false
	}
	screen = vec3.T{
		viewport.X + (ndcX+1)*0.5*viewport.Width,
		viewport.Y + (1-ndcY)*0.5*viewport.Height,
		(ndcZ + 1) * 0.5,
	}
	if !finiteFloat(screen[0]) || !finiteFloat(screen[1]) || !finiteFloat(screen[2]) {
		return vec3.Zero, false
	}
	return screen, true
}

// Unproject reverses a Project() call: it transforms the top-left origin
// pixel coordinates screen.X, screen.Y together with the [0, 1] depth
// screen.Z back into a world-space point using the given view and projection
// matrices and viewport.
//
// ok is true when screen corresponds to a point inside the clip volume
// (viewport rectangle with depth in [0, 1], with a small tolerance for
// rounding on the boundaries) and the combined view-projection matrix is
// invertible. ok is false for pixels outside of the viewport rectangle or
// depth range, for a singular view-projection matrix, for a nearly zero
// homogeneous w during the reverse transformation, for a zero or negative
// viewport width or height, and for any non-finite input. On failure the
// returned vector is vec3.Zero.
//
// float32 and float64 implementations of this function have identical
// semantics, only their precision differs.
func Unproject(screen *vec3.T, view, projection *T, viewport Viewport) (world vec3.T, ok bool) {
	if !finiteFloat(screen[0]) || !finiteFloat(screen[1]) || !finiteFloat(screen[2]) ||
		!finiteMat(view) || !finiteMat(projection) ||
		!finiteFloat(viewport.X) || !finiteFloat(viewport.Y) ||
		!finiteFloat(viewport.Width) || !finiteFloat(viewport.Height) {
		return vec3.Zero, false
	}
	if viewport.Width <= 0 || viewport.Height <= 0 {
		return vec3.Zero, false
	}

	var viewProjection T
	viewProjection.AssignMul(projection, view)
	inverse, invertible := invert4x4(&viewProjection)
	if !invertible {
		return vec3.Zero, false
	}

	// Pixel coordinates back to normalized device coordinates.
	ndcX := 2*(screen[0]-viewport.X)/viewport.Width - 1
	ndcY := 1 - 2*(screen[1]-viewport.Y)/viewport.Height
	ndcZ := 2*screen[2] - 1
	if !finiteFloat(ndcX) || !finiteFloat(ndcY) || !finiteFloat(ndcZ) ||
		outsideClipCube(ndcX, ndcY, ndcZ) {
		return vec3.Zero, false
	}

	// Reverse the perspective division of Project() by carrying w = 1 through
	// the inverse matrix and dividing by the resulting w.
	clip := inverse.MulVec4(&vec4.T{ndcX, ndcY, ndcZ, 1})
	if math.Abs(clip[3]) < minClipW {
		return vec3.Zero, false
	}
	world = vec3.T{
		clip[0] / clip[3],
		clip[1] / clip[3],
		clip[2] / clip[3],
	}
	if !finiteFloat(world[0]) || !finiteFloat(world[1]) || !finiteFloat(world[2]) {
		return vec3.Zero, false
	}
	return world, true
}
