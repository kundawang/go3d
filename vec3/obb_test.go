package vec3

import (
	stdmath "math"
	"testing"
)

const obbTestEps = 1e-4

// rotateZ45 rotates p by 45 degrees around the Z axis.
func rotateZ45(p *T) T {
	s := float32(stdmath.Sqrt(0.5))
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
			s := float32(1)
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
		if l := obb.Axis[a].Length(); stdmath.Abs(float64(l-1)) > obbTestEps {
			t.Errorf("axis %d not unit length: %v", a, l)
		}
		for b := a + 1; b < 3; b++ {
			if d := Dot(&obb.Axis[a], &obb.Axis[b]); stdmath.Abs(float64(d)) > obbTestEps {
				t.Errorf("axes %d and %d not orthogonal: %v", a, b, d)
			}
		}
	}

	// Axes must match the rotated canonical axes (up to sign and order).
	s := float32(stdmath.Sqrt(0.5))
	expected := [3]T{{s, s, 0}, {-s, s, 0}, {0, 0, 1}}
	for e := range expected {
		best := float32(0)
		for a := 0; a < 3; a++ {
			if d := stdmath.Abs(float64(Dot(&expected[e], &obb.Axis[a]))); d > float64(best) {
				best = float32(d)
			}
		}
		if stdmath.Abs(float64(best-1)) > obbTestEps {
			t.Errorf("expected axis %v not found, best |dot| = %v", expected[e], best)
		}
	}

	// Half sizes must be {2, 1, 0.5} as a set.
	wantHalf := []float32{2, 1, 0.5}
	gotHalf := []float32{obb.HalfSize[0], obb.HalfSize[1], obb.HalfSize[2]}
	for _, w := range wantHalf {
		found := false
		for i, g := range gotHalf {
			if stdmath.Abs(float64(g-w)) <= obbTestEps {
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
		if d := stdmath.Abs(float64(obb.Center[i] - offset[i])); d > obbTestEps {
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

func TestFitOBBTolerance(t *testing.T) {
	points := rotatedBoxCorners(T{1, 2, 3})
	const tol = 0.25
	obb := FitOBB(points, tol)
	// Every input point is contained even with zero query tolerance.
	for i := range points {
		if !obb.ContainsPoint(&points[i], 0) {
			t.Errorf("point %v not contained with fit tolerance", points[i])
		}
	}
	// Half sizes grew by exactly the tolerance.
	want := []float32{2 + tol, 1 + tol, 0.5 + tol}
	for _, w := range want {
		found := false
		for a := 0; a < 3; a++ {
			if stdmath.Abs(float64(obb.HalfSize[a]-w)) <= obbTestEps {
				found = true
				break
			}
		}
		if !found {
			t.Errorf("half size %v not found in %v", w, obb.HalfSize)
		}
	}
}

// obbEquivalent checks that two boxes describe the same volume: same center,
// same half sizes as a set, and axes equal up to sign and order.
func obbEquivalent(t *testing.T, a, b *OBB, eps float64, what string) {
	t.Helper()
	for i := 0; i < 3; i++ {
		if d := stdmath.Abs(float64(a.Center[i] - b.Center[i])); d > eps {
			t.Errorf("%s: center differs: %v vs %v", what, a.Center, b.Center)
			break
		}
	}
	used := [3]bool{}
	for i := 0; i < 3; i++ {
		best, bestJ := 0.0, -1
		for j := 0; j < 3; j++ {
			if used[j] {
				continue
			}
			if d := stdmath.Abs(float64(Dot(&a.Axis[i], &b.Axis[j]))); d > best {
				best, bestJ = d, j
			}
		}
		if bestJ < 0 || stdmath.Abs(best-1) > eps {
			t.Errorf("%s: axis %v has no match in %v", what, a.Axis[i], b.Axis)
			continue
		}
		used[bestJ] = true
		if d := stdmath.Abs(float64(a.HalfSize[i] - b.HalfSize[bestJ])); d > eps {
			t.Errorf("%s: half size %v vs %v along matched axes", what, a.HalfSize[i], b.HalfSize[bestJ])
		}
	}
}

func TestFitOBBInvariance(t *testing.T) {
	base := rotatedBoxCorners(T{0, 0, 0})
	// Add some interior points for a less symmetric point set.
	base = append(base, T{0.3, 0.1, 0.2}, T{-0.2, 0.4, -0.1})
	obbBase := FitOBB(base, 0)

	// Shuffled input.
	shuffled := make([]T, len(base))
	for i, p := range base {
		shuffled[(i*7+3)%len(base)] = p
	}
	obbShuffled := FitOBB(shuffled, 0)
	obbEquivalent(t, &obbBase, &obbShuffled, 1e-3, "shuffled")

	// Duplicated points.
	dup := append(append([]T{}, base...), base...)
	obbDup := FitOBB(dup, 0)
	obbEquivalent(t, &obbBase, &obbDup, 1e-3, "duplicated")

	// Translated input.
	offset := T{-7, 11, 2.5}
	translated := make([]T, len(base))
	for i, p := range base {
		translated[i] = *p.Add(&offset)
	}
	obbTranslated := FitOBB(translated, 0)
	wantCenter := obbBase.Center
	wantCenter.Add(&offset)
	for i := 0; i < 3; i++ {
		if d := stdmath.Abs(float64(obbTranslated.Center[i] - wantCenter[i])); d > 1e-3 {
			t.Errorf("translated: center %v, want %v", obbTranslated.Center, wantCenter)
			break
		}
	}
	obbTranslated.Center = obbBase.Center // compare the rest
	obbEquivalent(t, &obbBase, &obbTranslated, 1e-3, "translated")

	// Rotated input (30 degrees around Z): the fitted box must rotate along.
	const angle = stdmath.Pi / 6
	c, sn := float32(stdmath.Cos(angle)), float32(stdmath.Sin(angle))
	rot := func(p *T) T { return T{c*p[0] - sn*p[1], sn*p[0] + c*p[1], p[2]} }
	rotated := make([]T, len(base))
	for i, p := range base {
		rotated[i] = rot(&p)
	}
	obbRotated := FitOBB(rotated, 0)
	// Rotate the base fit the same way and compare.
	var want OBB
	want.Center = rot(&obbBase.Center)
	want.HalfSize = obbBase.HalfSize
	for a := 0; a < 3; a++ {
		want.Axis[a] = rot(&obbBase.Axis[a])
	}
	obbEquivalent(t, &want, &obbRotated, 1e-3, "rotated")
}

func TestFitOBBDegenerate(t *testing.T) {
	// Empty input.
	obb := FitOBB(nil, 0)
	if obb.Axis != [3]T{UnitX, UnitY, UnitZ} || !obb.Center.IsZero() || !obb.HalfSize.IsZero() {
		t.Errorf("empty input: got %+v", obb)
	}

	// Single point.
	p := T{1, 2, 3}
	obb = FitOBB([]T{p}, 0)
	if !obb.ContainsPoint(&p, 0) {
		t.Errorf("single point not contained")
	}

	// Collinear points along (1, 1, 1).
	dir := T{1, 1, 1}
	dir.Normalize()
	var line []T
	shift := T{4, -2, 7}
	for i := -5; i <= 5; i++ {
		p := dir.Scaled(float32(i))
		p.Add(&shift)
		line = append(line, p)
	}
	obb = FitOBB(line, 0)
	if d := stdmath.Abs(float64(Dot(&obb.Axis[0], &dir))); stdmath.Abs(d-1) > 1e-4 {
		t.Errorf("collinear: major axis %v not parallel to %v", obb.Axis[0], dir)
	}
	for a := 0; a < 3; a++ {
		if l := obb.Axis[a].Length(); stdmath.Abs(float64(l-1)) > 1e-4 {
			t.Errorf("collinear: axis %d not unit: %v", a, l)
		}
	}
	for i := range line {
		if !obb.ContainsPoint(&line[i], 1e-3) {
			t.Errorf("collinear: point %v not contained", line[i])
		}
	}
	// Deterministic under shuffling.
	shuffled := []T{line[3], line[9], line[0], line[7], line[1], line[10], line[5], line[2], line[8], line[4], line[6]}
	obb2 := FitOBB(shuffled, 0)
	for a := 0; a < 3; a++ {
		if d := Dot(&obb.Axis[a], &obb2.Axis[a]); stdmath.Abs(float64(d)-1) > 1e-3 {
			t.Errorf("collinear: axis %d not deterministic: %v vs %v", a, obb.Axis[a], obb2.Axis[a])
		}
	}

	// Coplanar points in the plane z = 5.
	var plane []T
	for x := -2; x <= 2; x++ {
		for y := -1; y <= 3; y++ {
			plane = append(plane, T{float32(x), float32(y), 5})
		}
	}
	obb = FitOBB(plane, 0)
	if d := stdmath.Abs(float64(obb.Axis[2][2])); stdmath.Abs(d-1) > 1e-4 {
		t.Errorf("coplanar: normal axis %v not parallel to Z", obb.Axis[2])
	}
	for i := range plane {
		if !obb.ContainsPoint(&plane[i], 1e-3) {
			t.Errorf("coplanar: point %v not contained", plane[i])
		}
	}

	// Near-equal major axes: square in the plane z = 0. The fit must be
	// deterministic under shuffling.
	square := []T{{1, 1, 0}, {-1, 1, 0}, {-1, -1, 0}, {1, -1, 0}}
	obbA := FitOBB(square, 0)
	obbB := FitOBB([]T{square[2], square[0], square[3], square[1]}, 0)
	for a := 0; a < 3; a++ {
		if d := Dot(&obbA.Axis[a], &obbB.Axis[a]); stdmath.Abs(float64(d)-1) > 1e-3 {
			t.Errorf("near-equal: axis %d not deterministic: %v vs %v", a, obbA.Axis[a], obbB.Axis[a])
		}
	}
	for i := range square {
		if !obbA.ContainsPoint(&square[i], 1e-3) {
			t.Errorf("near-equal: point %v not contained", square[i])
		}
	}
}
