package frustum

import (
	"errors"
	"math/rand"
	"testing"

	"math"
	"github.com/ungerik/go3d/float64/mat4"
	"github.com/ungerik/go3d/float64/vec3"
	"github.com/ungerik/go3d/float64/vec4"
)

// referenceClassify is an independent implementation used to cross-check
// the frustum classification: it transforms the 8 box corners into clip
// space and applies the homogeneous clip tests per corner. It also
// returns the minimum absolute plane distance encountered, so callers
// can skip borderline cases where different floating-point evaluation
// orders may legitimately disagree.
func referenceClassify(viewProj, transform *mat4.T, box *vec3.Box) (Relation, float64) {
	var m mat4.T
	m.AssignMul(viewProj, transform)

	minDist := float64(math.Inf(1))
	relation := Inside
	for plane := 0; plane < 6; plane++ {
		allOutside := true
		for i := 0; i < 8; i++ {
			corner := vec3.T{box.Min[0], box.Min[1], box.Min[2]}
			if i&1 != 0 {
				corner[0] = box.Max[0]
			}
			if i&2 != 0 {
				corner[1] = box.Max[1]
			}
			if i&4 != 0 {
				corner[2] = box.Max[2]
			}
			c4 := vec4.FromVec3(&corner)
			clip := m.MulVec4(&c4)
			dist := clipPlaneDistance(&clip, plane)
			if math.Abs(dist) < minDist {
				minDist = math.Abs(dist)
			}
			if dist >= 0 {
				allOutside = false
			} else {
				relation = Intersect
			}
		}
		if allOutside {
			return Outside, minDist
		}
	}
	return relation, minDist
}

// clipPlaneDistance returns the signed homogeneous distance of a clip
// space position to one clip plane: inside means -w <= x, y, z <= w.
func clipPlaneDistance(c *vec4.T, plane int) float64 {
	switch plane {
	case 0:
		return c[0] + c[3]
	case 1:
		return c[3] - c[0]
	case 2:
		return c[1] + c[3]
	case 3:
		return c[3] - c[1]
	case 4:
		return c[2] + c[3]
	default:
		return c[3] - c[2]
	}
}

func perspectiveViewProj(t *testing.T) *mat4.T {
	t.Helper()
	vp := new(mat4.T).AssignPerspective(math.Pi/3, 1, 1, 100)
	if _, err := FromViewProjection(vp); err != nil {
		t.Fatalf("perspective matrix rejected: %v", err)
	}
	return vp
}

func TestFromViewProjectionPerspective(t *testing.T) {
	f := MustFromViewProjection(perspectiveViewProj(t))
	// Planes must be normalized.
	for i, p := range f.Planes {
		l := math.Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
		if math.Abs(l-1) > 1e-5 {
			t.Errorf("plane %d not normalized: length %v", i, l)
		}
	}
	// A point in front of the camera must be inside, one behind outside.
	if !f.ContainsPoint(&vec3.T{0, 0, -10}) {
		t.Error("point in front of the camera should be inside the frustum")
	}
	if f.ContainsPoint(&vec3.T{0, 0, 10}) {
		t.Error("point behind the camera should be outside the frustum")
	}
}

func TestClassifyBoxPerspective(t *testing.T) {
	f := MustFromViewProjection(perspectiveViewProj(t))
	tests := []struct {
		name string
		box  vec3.Box
		want Relation
	}{
		{"centered in front", vec3.Box{Min: vec3.T{-1, -1, -10}, Max: vec3.T{1, 1, -5}}, Inside},
		{"behind camera", vec3.Box{Min: vec3.T{-1, -1, 1}, Max: vec3.T{1, 1, 5}}, Outside},
		{"beyond far plane", vec3.Box{Min: vec3.T{-1, -1, -200}, Max: vec3.T{1, 1, -150}}, Outside},
		{"left of frustum", vec3.Box{Min: vec3.T{-110, -1, -10}, Max: vec3.T{-100, 1, -5}}, Outside},
		{"right of frustum", vec3.Box{Min: vec3.T{100, -1, -10}, Max: vec3.T{110, 1, -5}}, Outside},
		{"above frustum", vec3.Box{Min: vec3.T{-1, 100, -10}, Max: vec3.T{1, 110, -5}}, Outside},
		{"straddling near plane", vec3.Box{Min: vec3.T{-0.5, -0.5, -2}, Max: vec3.T{0.5, 0.5, 0.5}}, Intersect},
		{"straddling far plane", vec3.Box{Min: vec3.T{-1, -1, -120}, Max: vec3.T{1, 1, -90}}, Intersect},
		// Half width at depth 10 with fovy 60 deg, aspect 1: tan(30)*10 ~ 5.77.
		{"crossing left plane", vec3.Box{Min: vec3.T{-8, -1, -11}, Max: vec3.T{-5, 1, -9}}, Intersect},
		// Zero-volume boxes behave like points.
		{"point inside", vec3.Box{Min: vec3.T{0, 0, -10}, Max: vec3.T{0, 0, -10}}, Inside},
		{"point outside", vec3.Box{Min: vec3.T{0, 0, 10}, Max: vec3.T{0, 0, 10}}, Outside},
		// Box with swapped Min/Max is ordered before classification.
		{"swapped min max", vec3.Box{Min: vec3.T{1, 1, -5}, Max: vec3.T{-1, -1, -10}}, Inside},
	}
	for _, tt := range tests {
		if got := f.ClassifyBox(&tt.box); got != tt.want {
			t.Errorf("%s: got %v, want %v", tt.name, got, tt.want)
		}
	}
}

// TestClassifyBoxTangent uses projection parameters that produce exact
// plane equations, so boxes exactly tangent to a plane can be tested:
// they must be classified Intersect, never Outside.
func TestClassifyBoxTangent(t *testing.T) {
	// Near plane z = -1 and far plane z = -3 have exact equations
	// (0,0,-1,-1) and (0,0,1,3) for znear=1, zfar=3.
	persp := new(mat4.T).AssignPerspective(math.Pi/3, 1, 1, 3)
	f := MustFromViewProjection(persp)
	tangentNear := vec3.Box{Min: vec3.T{-0.2, -0.2, -2}, Max: vec3.T{0.2, 0.2, -1}}
	if got := f.ClassifyBox(&tangentNear); got != Intersect {
		t.Errorf("box tangent to near plane: got %v, want Intersect", got)
	}
	tangentFar := vec3.Box{Min: vec3.T{-0.2, -0.2, -3}, Max: vec3.T{0.2, 0.2, -2}}
	if got := f.ClassifyBox(&tangentFar); got != Intersect {
		t.Errorf("box tangent to far plane: got %v, want Intersect", got)
	}

	// Orthogonal projection with exact plane equations:
	// right plane (-1,0,0,1), near plane (0,0,-1,-1).
	ortho := new(mat4.T).AssignOrthogonalProjection(-1, 1, -1, 1, 1, 3)
	fo := MustFromViewProjection(ortho)
	tangentRight := vec3.Box{Min: vec3.T{0.5, 0, -2}, Max: vec3.T{1, 0.5, -1.5}}
	if got := fo.ClassifyBox(&tangentRight); got != Intersect {
		t.Errorf("ortho box tangent to right plane: got %v, want Intersect", got)
	}
	tangentNearOrtho := vec3.Box{Min: vec3.T{0, 0, -2}, Max: vec3.T{0.5, 0.5, -1}}
	if got := fo.ClassifyBox(&tangentNearOrtho); got != Intersect {
		t.Errorf("ortho box tangent to near plane: got %v, want Intersect", got)
	}
	// Just outside the tangent box must still be culled.
	outsideRight := vec3.Box{Min: vec3.T{1.5, 0, -2}, Max: vec3.T{2, 0.5, -1.5}}
	if got := fo.ClassifyBox(&outsideRight); got != Outside {
		t.Errorf("ortho box right of right plane: got %v, want Outside", got)
	}
}

func TestClassifyBoxOrthogonal(t *testing.T) {
	vp := new(mat4.T).AssignOrthogonalProjection(-10, 10, -10, 10, 1, 100)
	f := MustFromViewProjection(vp)
	tests := []struct {
		name string
		box  vec3.Box
		want Relation
	}{
		{"centered", vec3.Box{Min: vec3.T{5, 5, -50}, Max: vec3.T{6, 6, -40}}, Inside},
		{"right of volume", vec3.Box{Min: vec3.T{20, 0, -50}, Max: vec3.T{30, 1, -40}}, Outside},
		{"behind camera", vec3.Box{Min: vec3.T{0, 0, 1}, Max: vec3.T{1, 1, 5}}, Outside},
		{"beyond far", vec3.Box{Min: vec3.T{0, 0, -200}, Max: vec3.T{1, 1, -150}}, Outside},
		{"crossing right edge", vec3.Box{Min: vec3.T{9, 0, -50}, Max: vec3.T{11, 1, -40}}, Intersect},
		{"crossing near plane", vec3.Box{Min: vec3.T{0, 0, -2}, Max: vec3.T{1, 1, 5}}, Intersect},
	}
	for _, tt := range tests {
		if got := f.ClassifyBox(&tt.box); got != tt.want {
			t.Errorf("%s: got %v, want %v", tt.name, got, tt.want)
		}
	}
}

// leftHandedPerspective returns a left-handed perspective projection
// (camera looks down +Z) mapping to OpenGL-style clip space.
func leftHandedPerspective(fovy, aspect, znear, zfar float64) mat4.T {
	tanHalf := math.Tan(fovy / 2)
	m := mat4.Zero
	m[0][0] = 1 / (aspect * tanHalf)
	m[1][1] = 1 / tanHalf
	m[2][2] = (zfar + znear) / (zfar - znear)
	m[3][2] = -2 * zfar * znear / (zfar - znear)
	m[2][3] = 1
	return m
}

func TestClassifyBoxLeftHanded(t *testing.T) {
	vp := leftHandedPerspective(math.Pi/3, 1, 1, 100)
	f := MustFromViewProjection(&vp)
	tests := []struct {
		name string
		box  vec3.Box
		want Relation
	}{
		{"in front (+Z)", vec3.Box{Min: vec3.T{-1, -1, 5}, Max: vec3.T{1, 1, 10}}, Inside},
		{"behind camera (-Z)", vec3.Box{Min: vec3.T{-1, -1, -10}, Max: vec3.T{1, 1, -5}}, Outside},
		{"beyond far", vec3.Box{Min: vec3.T{-1, -1, 150}, Max: vec3.T{1, 1, 200}}, Outside},
		{"straddling near", vec3.Box{Min: vec3.T{-0.5, -0.5, 0.5}, Max: vec3.T{0.5, 0.5, 2}}, Intersect},
	}
	for _, tt := range tests {
		if got := f.ClassifyBox(&tt.box); got != tt.want {
			t.Errorf("%s: got %v, want %v", tt.name, got, tt.want)
		}
	}
}

func TestClassifyTransformedBox(t *testing.T) {
	vp := perspectiveViewProj(t)
	f := MustFromViewProjection(vp)
	unit := vec3.Box{Min: vec3.T{-1, -1, -1}, Max: vec3.T{1, 1, 1}}

	tests := []struct {
		name      string
		box       vec3.Box
		transform *mat4.T
		want      Relation
	}{
		{
			"rotated 45 deg in front",
			unit,
			new(mat4.T).AssignYRotation(math.Pi/4).Translate(&vec3.T{0, 0, -10}),
			Inside,
		},
		{
			"rotated and behind camera",
			unit,
			new(mat4.T).AssignYRotation(math.Pi/4).Translate(&vec3.T{0, 0, 10}),
			Outside,
		},
		{
			// Rotating the unit box 45 deg around Y grows its footprint
			// to 2*sqrt(2) ~ 2.83; the rotated corners straddle the left
			// plane (half width ~5.77 at depth 10).
			"rotated corner crossing left plane",
			unit,
			new(mat4.T).AssignYRotation(math.Pi/4).Translate(&vec3.T{-6.2, 0, -10}),
			Intersect,
		},
		{
			"rotated and fully left of frustum",
			unit,
			new(mat4.T).AssignYRotation(math.Pi/4).Translate(&vec3.T{-20, 0, -10}),
			Outside,
		},
		{
			// Mirroring (negative scale) must not confuse the test.
			"mirrored in front",
			unit,
			new(mat4.T).SetScaling(&vec4.T{-1, 1, 1, 1}).Translate(&vec3.T{0, 0, -10}),
			Inside,
		},
		{
			// Singular transform: box collapsed onto the plane z = -10.
			"zero z scale inside",
			unit,
			new(mat4.T).SetScaling(&vec4.T{1, 1, 0, 1}).Translate(&vec3.T{0, 0, -10}),
			Inside,
		},
		{
			// Singular transform: box collapsed behind the camera.
			"zero z scale outside",
			unit,
			new(mat4.T).SetScaling(&vec4.T{1, 1, 0, 1}).Translate(&vec3.T{0, 0, 10}),
			Outside,
		},
		{
			// Zero-volume box transformed like a point.
			"degenerate point box in front",
			vec3.Box{Min: vec3.T{0, 0, 0}, Max: vec3.T{0, 0, 0}},
			new(mat4.T).AssignYRotation(math.Pi/4).Translate(&vec3.T{0, 0, -10}),
			Inside,
		},
	}
	for _, tt := range tests {
		if got := f.ClassifyTransformedBox(&tt.box, tt.transform); got != tt.want {
			t.Errorf("%s: got %v, want %v", tt.name, got, tt.want)
		}
		// Cross-check against the corner-based reference.
		if ref, _ := referenceClassify(vp, tt.transform, &tt.box); ref != tt.want {
			t.Errorf("%s: reference classification %v, want %v", tt.name, ref, tt.want)
		}
	}
}

func TestClassifyTransformedBoxRandomized(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	tested := 0
	for i := 0; i < 2000; i++ {
		var vp mat4.T
		if rng.Intn(2) == 0 {
			vp.AssignPerspective(0.3+rng.Float64()*2, 0.5+rng.Float64()*2, 0.1+rng.Float64()*10, 20+rng.Float64()*200)
		} else {
			vp.AssignOrthogonalProjection(-10, 10, -10, 10, 0.1+rng.Float64()*10, 20+rng.Float64()*200)
		}
		f := MustFromViewProjection(&vp)

		var box vec3.Box
		for c := 0; c < 3; c++ {
			a, b := rng.Float64()*4-2, rng.Float64()*4-2
			box.Min[c], box.Max[c] = math.Min(a, b), math.Max(a, b)
		}

		// Random affine transform: rotation, translation, and a scale
		// that may be negative (mirroring) or zero (singular).
		var transform mat4.T
		transform.AssignEulerRotation(rng.Float64()*math.Pi, rng.Float64()*math.Pi, rng.Float64()*math.Pi)
		scale := vec3.T{rng.Float64()*4 - 2, rng.Float64()*4 - 2, rng.Float64()*4 - 2}
		if i%7 == 0 {
			scale[rng.Intn(3)] = 0
		}
		transform.ScaleVec3(&scale)
		transform.Translate(&vec3.T{rng.Float64()*60 - 30, rng.Float64()*60 - 30, rng.Float64()*60 - 40})

		want, minDist := referenceClassify(&vp, &transform, &box)
		// Skip borderline cases where the different floating-point
		// evaluation orders of the two implementations may disagree.
		if minDist < 1e-3 {
			continue
		}
		if got := f.ClassifyTransformedBox(&box, &transform); got != want {
			t.Fatalf("iteration %d: got %v, want %v (box %v, transform %v)", i, got, want, box.String(), transform.String())
		}
		tested++
	}
	if tested < 1500 {
		t.Fatalf("too many borderline cases skipped: only %d of 2000 tested", tested)
	}
}

func TestDegenerateHandling(t *testing.T) {
	if _, err := FromViewProjection(&mat4.Zero); !errors.Is(err, ErrDegenerateMatrix) {
		t.Errorf("zero matrix: got %v, want ErrDegenerateMatrix", err)
	}

	nanMat := mat4.Ident
	nanMat[0][0] = math.NaN()
	if _, err := FromViewProjection(&nanMat); !errors.Is(err, ErrNonFiniteMatrix) {
		t.Errorf("NaN matrix: got %v, want ErrNonFiniteMatrix", err)
	}

	// znear == zfar produces infinities in the projection matrix.
	infMat := new(mat4.T).AssignPerspective(math.Pi/3, 1, 10, 10)
	if _, err := FromViewProjection(infMat); !errors.Is(err, ErrNonFiniteMatrix) {
		t.Errorf("znear == zfar: got %v, want ErrNonFiniteMatrix", err)
	}

	f := MustFromViewProjection(perspectiveViewProj(t))

	// Non-finite box coordinates are conservatively kept, never culled.
	nanBox := vec3.Box{Min: vec3.T{math.NaN(), 0, 0}, Max: vec3.T{1, 1, 1}}
	if got := f.ClassifyBox(&nanBox); got != Intersect {
		t.Errorf("NaN box: got %v, want Intersect", got)
	}
	infBox := vec3.Box{Min: vec3.T{0, 0, -10}, Max: vec3.T{float64(math.Inf(1)), 1, -5}}
	if got := f.ClassifyBox(&infBox); got != Intersect {
		t.Errorf("infinite box: got %v, want Intersect", got)
	}

	// Non-finite transform is conservatively kept, never culled.
	nanTransform := mat4.Ident
	nanTransform[3][0] = math.NaN()
	box := vec3.Box{Min: vec3.T{-1, -1, -1}, Max: vec3.T{1, 1, 1}}
	if got := f.ClassifyTransformedBox(&box, &nanTransform); got != Intersect {
		t.Errorf("NaN transform: got %v, want Intersect", got)
	}

	// Non-finite points are never contained.
	nanPoint := vec3.T{math.NaN(), 0, -10}
	if f.ContainsPoint(&nanPoint) {
		t.Error("NaN point should not be contained")
	}
}

func TestRelationString(t *testing.T) {
	if Outside.String() != "Outside" || Intersect.String() != "Intersect" || Inside.String() != "Inside" {
		t.Errorf("unexpected strings: %v %v %v", Outside, Intersect, Inside)
	}
}
