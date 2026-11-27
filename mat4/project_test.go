package mat4

import (
	stdmath "math"
	"testing"

	math "github.com/ungerik/go3d/fmath"
	"github.com/ungerik/go3d/vec3"
	"github.com/ungerik/go3d/vec4"
)

const (
	projectPixelEPS float32 = 1e-3
	projectWorldEPS float32 = 1e-3
)

func approxVec3(a, b vec3.T, eps float32) bool {
	return math.Abs(a[0]-b[0]) <= eps &&
		math.Abs(a[1]-b[1]) <= eps &&
		math.Abs(a[2]-b[2]) <= eps
}

func nan32() float32 { return float32(stdmath.NaN()) }
func inf32() float32 { return float32(stdmath.Inf(1)) }

// identity returns a fresh identity matrix, so tests do not depend on the
// package-level Ident staying unmodified.
func identity() T {
	return T{
		vec4.T{1, 0, 0, 0},
		vec4.T{0, 1, 0, 0},
		vec4.T{0, 0, 1, 0},
		vec4.T{0, 0, 0, 1},
	}
}

// referenceProject is an independent projection implementation written
// directly on top of the existing mat4/vec4 primitives, used to cross-check
// Project().
func referenceProject(world *vec3.T, view, projection *T, viewport Viewport) vec3.T {
	var viewProjection T
	viewProjection.AssignMul(projection, view)
	clip := viewProjection.MulVec4(&vec4.T{world[0], world[1], world[2], 1})
	ndc := clip.Vec3DividedByW()
	return vec3.T{
		viewport.X + (ndc[0]+1)*0.5*viewport.Width,
		viewport.Y + (1-ndc[1])*0.5*viewport.Height,
		(ndc[2] + 1) * 0.5,
	}
}

// TestProjectOrthographicCenter verifies the mandated case:
// the center world point of an orthographic camera must map to the viewport
// center with depth 0.5. For viewport (10, 20, 200, 100) the center pixel is
// (110, 70).
func TestProjectOrthographicCenter(t *testing.T) {
	view := identity()
	var projection T
	projection.AssignOrthogonalProjection(-50, 50, -25, 25, 1, 100)
	viewport := Viewport{X: 10, Y: 20, Width: 200, Height: 100}

	// go3d projections map the near plane to NDC z=-1 and the far plane to
	// NDC z=+1, hence (0, 0, -50.5) is the depth-0.5 midpoint in world space.
	center := vec3.T{0, 0, -50.5}
	screen, ok := Project(&center, &view, &projection, viewport)
	if !ok {
		t.Fatalf("Project returned ok=false for in-frustum center point")
	}
	want := vec3.T{110, 70, 0.5}
	if !approxVec3(screen, want, projectPixelEPS) {
		t.Errorf("ortho center: got %v, want %v", screen, want)
	}

	// Cross-check against the reference implementation on the raw matrices.
	ref := referenceProject(&center, &view, &projection, viewport)
	if !approxVec3(screen, ref, projectPixelEPS) {
		t.Errorf("ortho center differs from reference math: %v vs %v", screen, ref)
	}

	// Unprojecting the viewport center must recover the world point.
	back, ok := Unproject(&screen, &view, &projection, viewport)
	if !ok {
		t.Fatalf("Unproject returned ok=false for viewport center")
	}
	if !approxVec3(back, center, projectWorldEPS) {
		t.Errorf("ortho center round trip: got %v, want %v", back, center)
	}
}

// TestProjectOrthographicCorners maps the eight frustum corners and checks
// that the visible four land exactly on the viewport corners, with depths
// 0 and 1.
func TestProjectOrthographicCorners(t *testing.T) {
	view := identity()
	var projection T
	projection.AssignOrthogonalProjection(-50, 50, -25, 25, 1, 100)
	viewport := Viewport{X: 10, Y: 20, Width: 200, Height: 100}

	cases := []struct {
		world vec3.T
		want  vec3.T
	}{
		{vec3.T{50, 25, -1}, vec3.T{210, 20, 0}},     // near, top-right
		{vec3.T{-50, 25, -1}, vec3.T{10, 20, 0}},     // near, top-left
		{vec3.T{50, -25, -1}, vec3.T{210, 120, 0}},   // near, bottom-right
		{vec3.T{-50, -25, -1}, vec3.T{10, 120, 0}},   // near, bottom-left
		{vec3.T{50, 25, -100}, vec3.T{210, 20, 1}},   // far, top-right
		{vec3.T{-50, -25, -100}, vec3.T{10, 120, 1}}, // far, bottom-left
		{vec3.T{0, 0, -1}, vec3.T{110, 70, 0}},       // near center
		{vec3.T{0, 0, -100}, vec3.T{110, 70, 1}},     // far center
	}
	for i, c := range cases {
		screen, ok := Project(&c.world, &view, &projection, viewport)
		if !ok {
			t.Errorf("case %d: Project returned ok=false for frustum corner %v", i, c.world)
			continue
		}
		if !approxVec3(screen, c.want, projectPixelEPS) {
			t.Errorf("case %d: got %v, want %v", i, screen, c.want)
		}
		ref := referenceProject(&c.world, &view, &projection, viewport)
		if !approxVec3(screen, ref, projectPixelEPS) {
			t.Errorf("case %d: %v differs from reference math %v", i, screen, ref)
		}
		// Boundary pixels must unproject back to the original corner.
		back, ok := Unproject(&screen, &view, &projection, viewport)
		if !ok {
			t.Errorf("case %d: Unproject ok=false for boundary pixel %v", i, screen)
			continue
		}
		if !approxVec3(back, c.world, projectWorldEPS) {
			t.Errorf("case %d round trip: got %v, want %v", i, back, c.world)
		}
	}
}

// TestProjectUnprojectRoundTripPerspective exercises perspective projection,
// a non-square viewport and a rotated (rigid) view matrix over several
// in-frustum samples, in both projection directions.
func TestProjectUnprojectRoundTripPerspective(t *testing.T) {
	var view, projection T
	view.AssignYRotation(0.5)
	view.Translate(&vec3.T{2, -1, 0.5})
	projection.AssignPerspective(math.Pi/3, 2.0, 1, 100)
	viewport := Viewport{X: 10, Y: 20, Width: 200, Height: 100} // non-square, offset

	inverseView, invertible := invert4x4(&view)
	if !invertible {
		t.Fatal("rigid view matrix reported as singular")
	}
	// Camera-space points, all inside the frustum (camera looks down -Z,
	// near plane at z=-1, far plane at z=-100).
	cameraPoints := []vec3.T{
		{0, 0, -2},
		{1, -1, -5},
		{-3, 2, -15},
		{4, 3, -40},
		{-5, -4, -80},
		{0.5, 0.25, -99},
	}
	for i, camera := range cameraPoints {
		world := inverseView.MulVec3(&camera)

		screen, ok := Project(&world, &view, &projection, viewport)
		if !ok {
			t.Errorf("sample %d: Project ok=false for in-frustum point %v", i, world)
			continue
		}
		if screen[2] <= 0 || screen[2] >= 1 {
			t.Errorf("sample %d: depth %v outside (0, 1)", i, screen[2])
		}
		if screen[0] <= viewport.X || screen[0] >= viewport.X+viewport.Width ||
			screen[1] <= viewport.Y || screen[1] >= viewport.Y+viewport.Height {
			t.Errorf("sample %d: pixel %v outside viewport %+v", i, screen, viewport)
		}

		// A rotated view must produce the same pixels as projecting the
		// camera-space point with an identity view.
		ident := identity()
		screenCamera, okCamera := Project(&camera, &ident, &projection, viewport)
		if !okCamera || !approxVec3(screen, screenCamera, projectPixelEPS) {
			t.Errorf("sample %d: rotated view %v != identity view %v", i, screen, screenCamera)
		}
		ref := referenceProject(&world, &view, &projection, viewport)
		if !approxVec3(screen, ref, projectPixelEPS) {
			t.Errorf("sample %d: %v differs from reference math %v", i, screen, ref)
		}

		back, ok := Unproject(&screen, &view, &projection, viewport)
		if !ok {
			t.Errorf("sample %d: Unproject ok=false", i)
			continue
		}
		if !approxVec3(back, world, projectWorldEPS) {
			t.Errorf("sample %d round trip: got %v, want %v", i, back, world)
		}

		// And re-projecting the recovered point reproduces the pixel.
		reScreen, ok := Project(&back, &view, &projection, viewport)
		if !ok || !approxVec3(reScreen, screen, projectPixelEPS) {
			t.Errorf("sample %d re-project: got %v, want %v", i, reScreen, screen)
		}
	}
}

// TestProjectOrthographicRoundTrip covers orthographic projection with a
// rotated view and a non-square viewport. Samples are authored in camera
// space and transformed to world space with the inverse view, so they are
// guaranteed to lie inside the frustum.
func TestProjectOrthographicRoundTrip(t *testing.T) {
	var view, projection T
	view.AssignXRotation(-0.7)
	view.Translate(&vec3.T{-3, 1, 2})
	projection.AssignOrthogonalProjection(-20, 20, -10, 10, 0.5, 50)
	viewport := Viewport{X: 0, Y: 0, Width: 320, Height: 240}

	inverseView, invertible := invert4x4(&view)
	if !invertible {
		t.Fatal("rigid view matrix reported as singular")
	}
	cameraPoints := []vec3.T{
		{0, 0, -2},
		{10, 5, -4},
		{-15, -7, -20},
		{19.9, 9.9, -49.5},
	}
	for i, camera := range cameraPoints {
		world := inverseView.MulVec3(&camera)
		screen, ok := Project(&world, &view, &projection, viewport)
		if !ok {
			t.Errorf("sample %d: Project ok=false for %v", i, world)
			continue
		}
		back, ok := Unproject(&screen, &view, &projection, viewport)
		if !ok {
			t.Errorf("sample %d: Unproject ok=false", i)
			continue
		}
		if !approxVec3(back, world, projectWorldEPS) {
			t.Errorf("sample %d round trip: got %v, want %v", i, back, world)
		}
	}
}

// TestProjectOutsideClip checks points outside the clip volume.
func TestProjectOutsideClip(t *testing.T) {
	view := identity()
	var projection T
	projection.AssignPerspective(math.Pi/3, 1.0, 1, 100)
	viewport := Viewport{0, 0, 100, 100}

	outside := []vec3.T{
		{0, 0, -150}, // beyond the far plane
		{0, 0, -0.5}, // in front of the near plane
		{100, 0, -2}, // far outside the horizontal frustum
		{0, 100, -2}, // far outside the vertical frustum
	}
	for i, world := range outside {
		screen, ok := Project(&world, &view, &projection, viewport)
		if ok {
			t.Errorf("case %d: expected ok=false for outside point %v, got %v", i, world, screen)
		}
		if screen != vec3.Zero {
			t.Errorf("case %d: expected zero result on failure, got %v", i, screen)
		}
	}
}

// TestProjectNearZeroW verifies that a homogeneous clip-space w near zero
// (the camera plane of a perspective projection) yields a defined failure
// instead of a division by zero.
func TestProjectNearZeroW(t *testing.T) {
	view := identity()
	var projection T
	projection.AssignPerspective(math.Pi/3, 1.0, 1, 100)
	viewport := Viewport{0, 0, 100, 100}

	// The world origin sits exactly on the camera plane with an identity view,
	// so clip-space w is zero.
	origin := vec3.T{0, 0, 0}
	screen, ok := Project(&origin, &view, &projection, viewport)
	if ok {
		t.Errorf("expected ok=false for zero-w point, got %v", screen)
	}
	if screen != vec3.Zero {
		t.Errorf("expected zero result for zero-w point, got %v", screen)
	}

	// Points extremely close to the camera plane are rejected too.
	near := vec3.T{0, 0, -1e-9}
	if screen, ok := Project(&near, &view, &projection, viewport); ok || screen != vec3.Zero {
		t.Errorf("near-zero w: got (%v, %v), want (zero, false)", screen, ok)
	}
}

// TestUnprojectSingular checks that a singular view-projection matrix
// produces a defined failure.
func TestUnprojectSingular(t *testing.T) {
	viewport := Viewport{0, 0, 100, 100}
	screen := vec3.T{50, 50, 0.5}

	// The zero matrix is singular.
	if world, ok := Unproject(&screen, &Zero, &Zero, viewport); ok || world != vec3.Zero {
		t.Errorf("zero matrix: got (%v, %v), want (zero, false)", world, ok)
	}

	// A matrix with a zero projection row is rank deficient.
	var singular T
	singular[1] = vec4.T{}
	ident := identity()
	if world, ok := Unproject(&screen, &ident, &singular, viewport); ok || world != vec3.Zero {
		t.Errorf("rank-deficient matrix: got (%v, %v), want (zero, false)", world, ok)
	}
}

// TestInvert4x4 verifies the internal inversion helper against identity
// multiplication and known singular inputs.
func TestInvert4x4(t *testing.T) {
	ident := identity()
	inv, ok := invert4x4(&ident)
	if !ok {
		t.Fatal("identity matrix reported as singular")
	}
	if inv != ident {
		t.Errorf("inverse of identity: got %v", inv)
	}

	// Locally parsed invertible matrix (determinant 60.5, same values as the
	// fixture used by the other matrix tests, but parsed here so a mutation
	// of the package-level fixture cannot affect this test).
	original, err := Parse("3 1 0.5 0 2 5 2 0 1 6 7 0 2 100 1 1")
	if err != nil {
		t.Fatalf("parse test matrix: %v", err)
	}
	inv, ok = invert4x4(&original)
	if !ok {
		t.Fatal("invertible matrix reported as singular")
	}
	var product T
	product.AssignMul(&original, &inv)
	for col := 0; col < 4; col++ {
		for row := 0; row < 4; row++ {
			want := float32(0)
			if col == row {
				want = 1
			}
			if math.Abs(product[col][row]-want) > projectPixelEPS {
				t.Errorf("M * M^-1 not identity at [%d][%d]: %v", col, row, product[col][row])
			}
		}
	}

	if _, ok := invert4x4(&Zero); ok {
		t.Error("zero matrix reported as invertible")
	}
	var rankDeficient T
	rankDeficient[2] = vec4.T{}
	if _, ok := invert4x4(&rankDeficient); ok {
		t.Error("rank-deficient matrix reported as invertible")
	}
}

// TestProjectInvalidViewport covers zero-sized and negative viewports.
func TestProjectInvalidViewport(t *testing.T) {
	view := identity()
	var projection T
	projection.AssignPerspective(math.Pi/3, 1.0, 1, 100)
	world := vec3.T{0, 0, -5}
	screen := vec3.T{50, 50, 0.5}

	viewports := []Viewport{
		{0, 0, 0, 100},
		{0, 0, 100, 0},
		{0, 0, 0, 0},
		{10, 20, -200, 100},
		{10, 20, 200, -100},
	}
	for i, vp := range viewports {
		if got, ok := Project(&world, &view, &projection, vp); ok || got != vec3.Zero {
			t.Errorf("Project case %d: got (%v, %v), want (zero, false)", i, got, ok)
		}
		if got, ok := Unproject(&screen, &view, &projection, vp); ok || got != vec3.Zero {
			t.Errorf("Unproject case %d: got (%v, %v), want (zero, false)", i, got, ok)
		}
	}
}

// TestUnprojectOutsideViewport checks pixel and depth range validation.
func TestUnprojectOutsideViewport(t *testing.T) {
	view := identity()
	var projection T
	projection.AssignOrthogonalProjection(-50, 50, -25, 25, 1, 100)
	viewport := Viewport{X: 10, Y: 20, Width: 200, Height: 100}

	outside := []vec3.T{
		{9, 70, 0.5},    // left of X
		{211, 70, 0.5},  // right of X+Width
		{110, 19, 0.5},  // above Y
		{110, 121, 0.5}, // below Y+Height
		{110, 70, -0.1}, // depth below 0
		{110, 70, 1.1},  // depth above 1
	}
	for i, screen := range outside {
		if world, ok := Unproject(&screen, &view, &projection, viewport); ok || world != vec3.Zero {
			t.Errorf("case %d: got (%v, %v), want (zero, false)", i, world, ok)
		}
	}

	// Inclusive boundaries must succeed.
	for _, screen := range []vec3.T{
		{10, 20, 0},
		{210, 120, 1},
		{110, 70, 0.5},
	} {
		if _, ok := Unproject(&screen, &view, &projection, viewport); !ok {
			t.Errorf("boundary pixel %v rejected", screen)
		}
	}
}

// TestProjectNonFiniteInputs ensures NaN and infinity inputs have a defined
// result instead of propagating NaN.
func TestProjectNonFiniteInputs(t *testing.T) {
	view := identity()
	var projection T
	projection.AssignPerspective(math.Pi/3, 1.0, 1, 100)
	viewport := Viewport{0, 0, 100, 100}
	goodWorld := vec3.T{0, 0, -5}
	goodScreen := vec3.T{50, 50, 0.5}

	nanWorld := vec3.T{nan32(), 0, -5}
	infWorld := vec3.T{inf32(), 0, -5}
	if got, ok := Project(&nanWorld, &view, &projection, viewport); ok || got != vec3.Zero {
		t.Errorf("NaN world: got (%v, %v), want (zero, false)", got, ok)
	}
	if got, ok := Project(&infWorld, &view, &projection, viewport); ok || got != vec3.Zero {
		t.Errorf("Inf world: got (%v, %v), want (zero, false)", got, ok)
	}

	nanScreen := vec3.T{nan32(), 50, 0.5}
	infScreen := vec3.T{50, inf32(), 0.5}
	if got, ok := Unproject(&nanScreen, &view, &projection, viewport); ok || got != vec3.Zero {
		t.Errorf("NaN screen: got (%v, %v), want (zero, false)", got, ok)
	}
	if got, ok := Unproject(&infScreen, &view, &projection, viewport); ok || got != vec3.Zero {
		t.Errorf("Inf screen: got (%v, %v), want (zero, false)", got, ok)
	}

	// A non-finite matrix element invalidates the whole operation.
	badView := view
	badView[3][0] = nan32()
	if got, ok := Project(&goodWorld, &badView, &projection, viewport); ok || got != vec3.Zero {
		t.Errorf("NaN matrix in Project: got (%v, %v), want (zero, false)", got, ok)
	}
	if got, ok := Unproject(&goodScreen, &badView, &projection, viewport); ok || got != vec3.Zero {
		t.Errorf("NaN matrix in Unproject: got (%v, %v), want (zero, false)", got, ok)
	}

	nanViewport := Viewport{0, 0, nan32(), 100}
	if got, ok := Project(&goodWorld, &view, &projection, nanViewport); ok || got != vec3.Zero {
		t.Errorf("NaN viewport in Project: got (%v, %v), want (zero, false)", got, ok)
	}
	if got, ok := Unproject(&goodScreen, &view, &projection, nanViewport); ok || got != vec3.Zero {
		t.Errorf("NaN viewport in Unproject: got (%v, %v), want (zero, false)", got, ok)
	}
}
