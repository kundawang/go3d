package vec3

import (
	stdmath "math"
	"testing"

	math "github.com/ungerik/go3d/fmath"
)

const capsuleTestEpsilon = float32(1e-5)

func vecPracticallyEquals(t *testing.T, name string, got, want T, epsilon float32) {
	t.Helper()
	for i := 0; i < 3; i++ {
		if math.Abs(got[i]-want[i]) > epsilon {
			t.Errorf("%s[%d]: got %v, want %v", name, i, got, want)
			return
		}
	}
}

func scalarPracticallyEquals(t *testing.T, name string, got, want, epsilon float32) {
	t.Helper()
	if math.Abs(got-want) > epsilon {
		t.Errorf("%s: got %v, want %v", name, got, want)
	}
}

// Unit box, capsule whose axis runs parallel to a face at distance 0.3,
// radius 0.4 (penetration 0.1), plus a smaller-radius non-contact control.
func TestCapsuleFace(t *testing.T) {
	box := Box{T{0, 0, 0}, T{1, 1, 1}}
	a := T{0.2, -0.3, 0.5}
	b := T{0.8, -0.3, 0.5}

	hit := box.IntersectCapsule(&a, &b, 0.4)
	if !hit.Contact {
		t.Fatal("expected contact")
	}
	scalarPracticallyEquals(t, "depth", hit.Depth, 0.1, capsuleTestEpsilon)
	vecPracticallyEquals(t, "normal", hit.Normal, T{0, -1, 0}, capsuleTestEpsilon)
	vecPracticallyEquals(t, "box point", hit.BoxPoint, T{0.5, 0, 0.5}, capsuleTestEpsilon)
	// Closest capsule surface point: axis y=-0.3 plus the radius toward the
	// box (normal -y, so the surface top is at y=-0.3+0.4=+0.1).
	vecPracticallyEquals(t, "capsule point", hit.CapsulePoint, T{0.5, 0.1, 0.5}, capsuleTestEpsilon)
	scalarPracticallyEquals(t, "normal length", hit.Normal.Length(), 1, capsuleTestEpsilon)

	// Smaller radius control: gap 0.1, no contact but full closest info.
	near := box.IntersectCapsule(&a, &b, 0.2)
	if near.Contact {
		t.Fatal("expected no contact for r=0.2")
	}
	if near.Depth != 0 {
		t.Errorf("non-contact depth must be zero, got %v", near.Depth)
	}
	vecPracticallyEquals(t, "near normal", near.Normal, T{0, -1, 0}, capsuleTestEpsilon)
	vecPracticallyEquals(t, "near box point", near.BoxPoint, T{0.5, 0, 0.5}, capsuleTestEpsilon)
	// Surface top toward the box: y=-0.3+0.2=-0.1.
	vecPracticallyEquals(t, "near capsule point", near.CapsulePoint, T{0.5, -0.1, 0.5}, capsuleTestEpsilon)

	// Tangent: distance == radius counts as contact with zero depth.
	tangent := box.IntersectCapsule(&a, &b, 0.3)
	if !tangent.Contact {
		t.Fatal("tangent configuration must count as contact")
	}
	scalarPracticallyEquals(t, "tangent depth", tangent.Depth, 0, capsuleTestEpsilon)
}

// Capsule facing a box corner diagonally: end-cap and box corner geometry.
func TestCapsuleCorner(t *testing.T) {
	box := Box{T{0, 0, 0}, T{1, 1, 1}}
	corner := T{1, 1, 1}
	a := T{1.3, 1.3, 1.3}
	b := T{3, 3, 3}
	inv := float32(1) / float32(math.Sqrt(3))
	diag := T{inv, inv, inv}

	// r=0.4: gap 0.3*sqrt(3) ≈ 0.5196 > r, no contact.
	hit := box.IntersectCapsule(&a, &b, 0.4)
	if hit.Contact {
		t.Fatal("expected no contact")
	}
	vecPracticallyEquals(t, "box point", hit.BoxPoint, corner, capsuleTestEpsilon)
	vecPracticallyEquals(t, "normal", hit.Normal, diag, capsuleTestEpsilon)
	// Capsule surface point is the axis point shifted by r toward the box.
	wantNear := T{1.3 - 0.4*inv, 1.3 - 0.4*inv, 1.3 - 0.4*inv}
	vecPracticallyEquals(t, "capsule point", hit.CapsulePoint, wantNear, capsuleTestEpsilon)
	scalarPracticallyEquals(t, "surface separation", Distance(&hit.BoxPoint, &hit.CapsulePoint),
		float32(math.Sqrt(3))*0.3-0.4, capsuleTestEpsilon)

	// r=0.6: same geometry now penetrates through the end cap.
	hit = box.IntersectCapsule(&a, &b, 0.6)
	if !hit.Contact {
		t.Fatal("expected contact")
	}
	scalarPracticallyEquals(t, "depth", hit.Depth, 0.6-float32(math.Sqrt(3))*0.3, capsuleTestEpsilon)
	vecPracticallyEquals(t, "normal", hit.Normal, diag, capsuleTestEpsilon)
	vecPracticallyEquals(t, "box point", hit.BoxPoint, corner, capsuleTestEpsilon)
	scalarPracticallyEquals(t, "normal length", hit.Normal.Length(), 1, capsuleTestEpsilon)

	// Smaller-radius control at distance 0.1 from the corner.
	c := T{1.1, 1.1, 1.1}
	near := box.IntersectCapsule(&c, &c, 0.05)
	if near.Contact {
		t.Fatal("expected no contact")
	}
	scalarPracticallyEquals(t, "near surface gap", Distance(&near.BoxPoint, &near.CapsulePoint),
		float32(math.Sqrt(3))*0.1-0.05, capsuleTestEpsilon)
}

// End cap of a capsule resting on a box edge.
func TestCapsuleEndCapEdge(t *testing.T) {
	box := Box{T{0, 0, 0}, T{1, 1, 1}}
	a := T{0.5, -0.25, -0.25}
	b := T{0.5, 1.5, -0.25}

	// Closest features: every segment point with y in [0,1] is at the same
	// 0.25 z-gap, so the flat minimum region is centered on t=0.5 giving the
	// segment point (0.5,0.5,-0.25) and the box edge point (0.5,0.5,0).
	// Radius 0.4 gives a 0.15 penetration.
	hit := box.IntersectCapsule(&a, &b, 0.4)
	if !hit.Contact {
		t.Fatal("expected contact")
	}
	scalarPracticallyEquals(t, "depth", hit.Depth, 0.15, capsuleTestEpsilon)
	vecPracticallyEquals(t, "normal", hit.Normal, T{0, 0, -1}, capsuleTestEpsilon)
	vecPracticallyEquals(t, "box point", hit.BoxPoint, T{0.5, 0.5, 0}, capsuleTestEpsilon)
	// Axis point z=-0.25 plus radius toward the box (normal -z): z=+0.15.
	vecPracticallyEquals(t, "capsule point", hit.CapsulePoint, T{0.5, 0.5, 0.15}, capsuleTestEpsilon)
}

// Capsule axis passing through the box volume: no unique closest feature,
// the capsule is pushed out through the face with the smallest overlap.
func TestCapsuleAxisThroughBox(t *testing.T) {
	box := Box{T{0, 0, 0}, T{1, 1, 1}}
	// Segment starts inside the box and exits through the +x face.
	a := T{0.7, 0.5, 0.5}
	b := T{2.0, 0.5, 0.5}

	hit := box.IntersectCapsule(&a, &b, 0.2)
	if !hit.Contact {
		t.Fatal("expected contact")
	}
	// Capsule x extent is [0.5,2.2]; only the +x face is crossed, overlap
	// 2.2-1.0 = 1.2. The y/z extents stay inside the box.
	scalarPracticallyEquals(t, "depth", hit.Depth, 1.2, capsuleTestEpsilon)
	vecPracticallyEquals(t, "normal", hit.Normal, T{1, 0, 0}, capsuleTestEpsilon)
	vecPracticallyEquals(t, "box point", hit.BoxPoint, T{1, 0.5, 0.5}, capsuleTestEpsilon)
	// Responsible endpoint b at x=2 plus radius along +x => x=2.2.
	vecPracticallyEquals(t, "capsule point", hit.CapsulePoint, T{2.2, 0.5, 0.5}, capsuleTestEpsilon)

	// Segment fully spanning the box along y; the two y faces overlap by
	// the same amount and the deterministic tie break picks the low face.
	c := T{0.5, -0.5, 0.5}
	d := T{0.5, 1.5, 0.5}
	hit = box.IntersectCapsule(&c, &d, 0.1)
	if !hit.Contact {
		t.Fatal("expected contact")
	}
	// Capsule y extent is [-0.6,1.6], overlap 0.6 on both y faces.
	scalarPracticallyEquals(t, "sym depth", hit.Depth, 0.6, capsuleTestEpsilon)
	vecPracticallyEquals(t, "sym normal", hit.Normal, T{0, -1, 0}, capsuleTestEpsilon)
	vecPracticallyEquals(t, "sym box point", hit.BoxPoint, T{0.5, 0, 0.5}, capsuleTestEpsilon)
	vecPracticallyEquals(t, "sym capsule point", hit.CapsulePoint, T{0.5, -0.6, 0.5}, capsuleTestEpsilon)
}

// A zero-length axis degenerates the capsule to a sphere.
func TestCapsuleZeroLengthIsSphere(t *testing.T) {
	box := Box{T{0, 0, 0}, T{1, 1, 1}}

	// Sphere center outside the box, face case.
	center := T{0.5, 0.7, 0.5}
	hit := box.IntersectCapsule(&center, &center, 0.4)
	if !hit.Contact {
		t.Fatal("expected contact")
	}
	scalarPracticallyEquals(t, "depth", hit.Depth, 0.1, capsuleTestEpsilon)
	vecPracticallyEquals(t, "normal", hit.Normal, T{0, 1, 0}, capsuleTestEpsilon)
	vecPracticallyEquals(t, "box point", hit.BoxPoint, T{0.5, 1, 0.5}, capsuleTestEpsilon)

	// Sphere center inside the box: push out through the nearest face.
	// Sphere x extent is [-0.2,0.6], low-x overlap is 0.2, all other
	// extents stay inside the box.
	inside := T{0.2, 0.5, 0.5}
	hit = box.IntersectCapsule(&inside, &inside, 0.4)
	if !hit.Contact {
		t.Fatal("expected contact")
	}
	scalarPracticallyEquals(t, "inside depth", hit.Depth, 0.2, capsuleTestEpsilon)
	vecPracticallyEquals(t, "inside normal", hit.Normal, T{-1, 0, 0}, capsuleTestEpsilon)
	vecPracticallyEquals(t, "inside box point", hit.BoxPoint, T{0, 0.5, 0.5}, capsuleTestEpsilon)
	vecPracticallyEquals(t, "inside capsule point", hit.CapsulePoint, T{-0.2, 0.5, 0.5}, capsuleTestEpsilon)
}

// Invalid and degenerate inputs must return a stable zero result.
func TestCapsuleInvalidInput(t *testing.T) {
	good := T{0, 0, 0}
	other := T{1, 1, 1}

	badBox := Box{T{1, 0, 0}, T{0, 1, 1}}
	if hit := badBox.IntersectCapsule(&good, &other, 0.1); hit.Contact {
		t.Error("invalid box must not report contact")
	}

	box := Box{good, other}
	if hit := box.IntersectCapsule(&good, &other, -0.1); hit.Contact {
		t.Error("negative radius must not report contact")
	}

	nan := T{float32(stdmath.NaN()), 0, 0}
	if hit := box.IntersectCapsule(&nan, &other, 0.1); hit.Contact {
		t.Error("NaN input must not report contact")
	}
	inf := T{float32(stdmath.Inf(1)), 0, 0}
	if hit := box.IntersectCapsule(&inf, &other, 0.1); hit.Contact {
		t.Error("Inf input must not report contact")
	}

	// Zero radius against a point on the surface: tangent contact.
	surface := T{0.5, 1, 0.5}
	hit := box.IntersectCapsule(&surface, &surface, 0)
	if !hit.Contact || hit.Depth != 0 {
		t.Errorf("zero-radius tangent must report zero-depth contact, got %+v", hit)
	}
}

// A segment grazing an edge tangentially must report contact with zero
// depth and a deterministic edge normal (smallest axis index wins).
func TestCapsuleTangentEdge(t *testing.T) {
	box := Box{T{0, 0, 0}, T{1, 1, 1}}
	// Parallel to the x-axis edge at y=0, z=0, offset in y and z equally.
	a := T{-1, -0.2, -0.2}
	b := T{2, -0.2, -0.2}
	hit := box.IntersectCapsule(&a, &b, float32(math.Sqrt(0.08)))
	if !hit.Contact {
		t.Fatal("tangent edge configuration must count as contact")
	}
	scalarPracticallyEquals(t, "depth", hit.Depth, 0, capsuleTestEpsilon)
	inv := float32(1) / float32(math.Sqrt(2))
	vecPracticallyEquals(t, "normal", hit.Normal, T{0, -inv, -inv}, capsuleTestEpsilon)
}
