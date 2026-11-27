package vec3

import (
	"math"
	"testing"
)

const obbTestEps = 1e-9

// rotateZ45 rotates p by 45 degrees around the Z axis.
func rotateZ45(p *T) T {
	s := math.Sqrt(0.5)
	return T{(p[0] - p[1]) * s, (p[0] + p[1]) * s, p[2]}
}

// rotatedBoxCorners returns the 8 corners of a 4x2x1 box rotated 45 degrees
// around Z and translated by offset.
func rotatedBoxCorners(offset T) []T {
	half := T{2, 1, 0.5}
	points := make([]T, 0, 8)
	for i := 0; i < 8; i++ {
		var p T
		for a := 0; a < 3; a++ {
			s := float64(1)
			if i&(1<<uint(a)) != 0 {
				s = -1
			}
			p[a] = s * half[a]
		}
		r := rotateZ45(&p)
		r.Add(&offset)
		points = append(points, r)
	}
	return points
}

func TestFitOBBRotatedBox(t *testing.T) {
	offset := T{10, -5, 3}
	points := rotatedBoxCorners(offset)
	obb := FitOBB(points, 0)

	// Axes must be unit length and pairwise orthogonal.
	for a := 0; a < 3; a++ {
		if l := obb.Axis[a].Length(); math.Abs(l-1) > obbTestEps {
			t.Errorf("axis %d not unit length: %v", a, l)
		}
		for b := a + 1; b < 3; b++ {
			if d := Dot(&obb.Axis[a], &obb.Axis[b]); math.Abs(d) > obbTestEps {
				t.Errorf("axes %d and %d not orthogonal: %v", a, b, d)
			}
		}
	}

	// Axes must match the rotated canonical axes (up to sign and order).
	s := math.Sqrt(0.5)
	expected := [3]T{{s, s, 0}, {-s, s, 0}, {0, 0, 1}}
	for e := range expected {
		best := 0.0
		for a := 0; a < 3; a++ {
			if d := math.Abs(Dot(&expected[e], &obb.Axis[a])); d > best {
				best = d
			}
		}
		if math.Abs(best-1) > obbTestEps {
			t.Errorf("expected axis %v not found, best |dot| = %v", expected[e], best)
		}
	}

	// Half sizes must be {2, 1, 0.5} as a set.
	wantHalf := []float64{2, 1, 0.5}
	gotHalf := []float64{obb.HalfSize[0], obb.HalfSize[1], obb.HalfSize[2]}
	for _, w := range wantHalf {
		found := false
		for i, g := range gotHalf {
			if math.Abs(g-w) <= obbTestEps {
				gotHalf = append(gotHalf[:i], gotHalf[i+1:]...)
				found = true
				break
			}
		}
		if !found {
			t.Errorf("half size %v not found in %v", w, obb.HalfSize)
		}
	}

	// Center must be the rotated box center, which is the offset.
	for i := 0; i < 3; i++ {
		if d := math.Abs(obb.Center[i] - offset[i]); d > obbTestEps {
			t.Errorf("center[%d] = %v, want %v", i, obb.Center[i], offset[i])
		}
	}

	// All input points must be contained.
	for i := range points {
		if !obb.ContainsPoint(&points[i], obbTestEps) {
			t.Errorf("point %v not contained", points[i])
		}
	}

	// A clearly outside point must not be contained.
	outside := obb.Center
	off := obb.Axis[0].Scaled(obb.HalfSize[0] * 3)
	outside.Add(&off)
	if obb.ContainsPoint(&outside, obbTestEps) {
		t.Errorf("outside point %v reported as contained", outside)
	}

	// The reconstructed corners must match the input corners as a set.
	corners := obb.Corners()
	for i := range points {
		found := false
		for j := range corners {
			if Distance(&points[i], &corners[j]) <= obbTestEps {
				found = true
				break
			}
		}
		if !found {
			t.Errorf("input corner %v not reproduced by Corners()", points[i])
		}
	}
}

func TestFitOBBDegenerate(t *testing.T) {
	// Collinear points along (1, 1, 1).
	dir := T{1, 1, 1}
	dir.Normalize()
	var line []T
	shift := T{4, -2, 7}
	for i := -5; i <= 5; i++ {
		p := dir.Scaled(float64(i))
		p.Add(&shift)
		line = append(line, p)
	}
	obb := FitOBB(line, 0)
	if d := math.Abs(Dot(&obb.Axis[0], &dir)); math.Abs(d-1) > 1e-9 {
		t.Errorf("collinear: major axis %v not parallel to %v", obb.Axis[0], dir)
	}
	for i := range line {
		if !obb.ContainsPoint(&line[i], 1e-9) {
			t.Errorf("collinear: point %v not contained", line[i])
		}
	}

	// Coplanar points in the plane z = 5.
	var plane []T
	for x := -2; x <= 2; x++ {
		for y := -1; y <= 3; y++ {
			plane = append(plane, T{float64(x), float64(y), 5})
		}
	}
	obb = FitOBB(plane, 0)
	if d := math.Abs(obb.Axis[2][2]); math.Abs(d-1) > 1e-9 {
		t.Errorf("coplanar: normal axis %v not parallel to Z", obb.Axis[2])
	}
	for i := range plane {
		if !obb.ContainsPoint(&plane[i], 1e-9) {
			t.Errorf("coplanar: point %v not contained", plane[i])
		}
	}
}
