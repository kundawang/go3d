package vec3

import (
	"errors"
	"math"
	"testing"
)

// sphereTestCubeCorners returns the 8 corners of the axis aligned cube
// [-halfExtent, halfExtent]^3 translated by offset.
func sphereTestCubeCorners(halfExtent float64, offset T) []T {
	corners := make([]T, 0, 8)
	for _, x := range []float64{-halfExtent, halfExtent} {
		for _, y := range []float64{-halfExtent, halfExtent} {
			for _, z := range []float64{-halfExtent, halfExtent} {
				corners = append(corners, Add(&T{x, y, z}, &offset))
			}
		}
	}
	return corners
}

func assertSphereContainsAllPoints(t *testing.T, s *Sphere, points []T) {
	t.Helper()
	for i := range points {
		if !s.ContainsPoint(&points[i]) {
			t.Errorf("point %v not contained in sphere %v", points[i], s)
		}
	}
}

func TestMinimumBoundingSphereCube(t *testing.T) {
	corners := sphereTestCubeCorners(1, T{})
	s, err := MinimumBoundingSphere(corners)
	if err != nil {
		t.Fatal(err)
	}
	if !s.Center.PracticallyEquals(&Zero, 1e-9) {
		t.Errorf("expected center %v, got %v", Zero, s.Center)
	}
	if expected := math.Sqrt(3); math.Abs(s.Radius-expected) > 1e-9 {
		t.Errorf("expected radius %v, got %v", expected, s.Radius)
	}
	assertSphereContainsAllPoints(t, &s, corners)
}

func TestMinimumBoundingSphereTranslatedCube(t *testing.T) {
	offset := T{10, -20, 30}
	corners := sphereTestCubeCorners(2, offset)
	s, err := MinimumBoundingSphere(corners)
	if err != nil {
		t.Fatal(err)
	}
	if !s.Center.PracticallyEquals(&offset, 1e-9) {
		t.Errorf("expected center %v, got %v", offset, s.Center)
	}
	if expected := 2 * math.Sqrt(3); math.Abs(s.Radius-expected) > 1e-9 {
		t.Errorf("expected radius %v, got %v", expected, s.Radius)
	}
	assertSphereContainsAllPoints(t, &s, corners)
}

func TestMinimumBoundingSphereCubePointOrder(t *testing.T) {
	corners := sphereTestCubeCorners(1, T{1, 2, 3})
	reference, err := MinimumBoundingSphere(corners)
	if err != nil {
		t.Fatal(err)
	}

	reversed := make([]T, len(corners))
	for i, p := range corners {
		reversed[len(corners)-1-i] = p
	}
	permuted := make([]T, len(corners))
	for i, p := range corners {
		permuted[(i*5+3)%len(corners)] = p
	}

	for name, reordered := range map[string][]T{"reversed": reversed, "permuted": permuted} {
		s, err := MinimumBoundingSphere(reordered)
		if err != nil {
			t.Fatal(err)
		}
		if !s.Center.PracticallyEquals(&reference.Center, 1e-9) {
			t.Errorf("%s order: expected center %v, got %v", name, reference.Center, s.Center)
		}
		if math.Abs(s.Radius-reference.Radius) > 1e-9 {
			t.Errorf("%s order: expected radius %v, got %v", name, reference.Radius, s.Radius)
		}
		assertSphereContainsAllPoints(t, &s, reordered)
	}
}

func TestMinimumBoundingSphereDuplicatePoints(t *testing.T) {
	corners := sphereTestCubeCorners(1, T{})
	duplicated := make([]T, 0, 3*len(corners))
	for i := 0; i < 3; i++ {
		duplicated = append(duplicated, corners...)
	}
	s, err := MinimumBoundingSphere(duplicated)
	if err != nil {
		t.Fatal(err)
	}
	if !s.Center.PracticallyEquals(&Zero, 1e-9) {
		t.Errorf("expected center %v, got %v", Zero, s.Center)
	}
	if expected := math.Sqrt(3); math.Abs(s.Radius-expected) > 1e-9 {
		t.Errorf("expected radius %v, got %v", expected, s.Radius)
	}
	assertSphereContainsAllPoints(t, &s, duplicated)
}

func TestMinimumBoundingSphereEmptyPointSet(t *testing.T) {
	if _, err := MinimumBoundingSphere(nil); !errors.Is(err, ErrEmptyPointSet) {
		t.Errorf("expected ErrEmptyPointSet, got %v", err)
	}
	if _, err := MinimumBoundingSphere([]T{}); !errors.Is(err, ErrEmptyPointSet) {
		t.Errorf("expected ErrEmptyPointSet, got %v", err)
	}
}

func TestMinimumBoundingSphereSinglePoint(t *testing.T) {
	point := T{1, 2, 3}
	s, err := MinimumBoundingSphere([]T{point})
	if err != nil {
		t.Fatal(err)
	}
	if s.Center != point {
		t.Errorf("expected center %v, got %v", point, s.Center)
	}
	if s.Radius != 0 {
		t.Errorf("expected radius 0, got %v", s.Radius)
	}
}

func TestMinimumBoundingSphereCollinearPoints(t *testing.T) {
	points := []T{{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}, {-4, 0, 0}}
	s, err := MinimumBoundingSphere(points)
	if err != nil {
		t.Fatal(err)
	}
	expectedCenter := T{-0.5, 0, 0}
	if !s.Center.PracticallyEquals(&expectedCenter, 1e-9) {
		t.Errorf("expected center %v, got %v", expectedCenter, s.Center)
	}
	if math.Abs(s.Radius-3.5) > 1e-9 {
		t.Errorf("expected radius 3.5, got %v", s.Radius)
	}
	assertSphereContainsAllPoints(t, &s, points)
}

func TestMinimumBoundingSphereCoplanarPoints(t *testing.T) {
	// Square in the plane z = 2.
	points := []T{{-1, -1, 2}, {1, -1, 2}, {1, 1, 2}, {-1, 1, 2}}
	s, err := MinimumBoundingSphere(points)
	if err != nil {
		t.Fatal(err)
	}
	expectedCenter := T{0, 0, 2}
	if !s.Center.PracticallyEquals(&expectedCenter, 1e-9) {
		t.Errorf("expected center %v, got %v", expectedCenter, s.Center)
	}
	if expected := math.Sqrt(2); math.Abs(s.Radius-expected) > 1e-9 {
		t.Errorf("expected radius %v, got %v", expected, s.Radius)
	}
	assertSphereContainsAllPoints(t, &s, points)
}

func TestMinimumBoundingSphereNearDegeneratePoints(t *testing.T) {
	// Square with one corner lifted by a tiny amount: almost coplanar.
	points := []T{{-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}, {0, 0, 1e-13}}
	s, err := MinimumBoundingSphere(points)
	if err != nil {
		t.Fatal(err)
	}
	if s.Radius != s.Radius { // NaN check
		t.Fatalf("radius is not a number: %v", s.Radius)
	}
	if expected := math.Sqrt(2); math.Abs(s.Radius-expected) > 0.05*expected {
		t.Errorf("expected radius close to %v, got %v", expected, s.Radius)
	}
	assertSphereContainsAllPoints(t, &s, points)
}

func TestMinimumBoundingSphereDoesNotModifyInput(t *testing.T) {
	points := sphereTestCubeCorners(1, T{3, 2, 1})
	backup := make([]T, len(points))
	copy(backup, points)
	if _, err := MinimumBoundingSphere(points); err != nil {
		t.Fatal(err)
	}
	for i := range points {
		if points[i] != backup[i] {
			t.Fatalf("input point %d modified: %v != %v", i, points[i], backup[i])
		}
	}
}

func TestMinimumBoundingSpherePointCloud(t *testing.T) {
	// Deterministic pseudo random point cloud in [-1, 1]^3.
	var state uint64 = 12345
	next := func() float64 {
		state = state*6364136223846793005 + 1442695040888963407
		return float64(state>>11)/float64(1<<53)*2 - 1
	}
	points := make([]T, 200)
	for i := range points {
		points[i] = T{next(), next(), next()}
	}

	s, err := MinimumBoundingSphere(points)
	if err != nil {
		t.Fatal(err)
	}
	assertSphereContainsAllPoints(t, &s, points)

	// The radius cannot be smaller than half the maximum pairwise distance.
	maxPairDistance := float64(0)
	for i := range points {
		for j := i + 1; j < len(points); j++ {
			if d := Distance(&points[i], &points[j]); d > maxPairDistance {
				maxPairDistance = d
			}
		}
	}
	if s.Radius < maxPairDistance/2-1e-9 {
		t.Errorf("radius %v smaller than half the maximum pairwise distance %v", s.Radius, maxPairDistance/2)
	}

	// The radius cannot exceed half the diagonal of the points' Box.
	box := Box{Min: points[0], Max: points[0]}
	for i := range points {
		box.Min = Min(&box.Min, &points[i])
		box.Max = Max(&box.Max, &points[i])
	}
	diagonal := box.Diagonal()
	if maxRadius := diagonal.Length() / 2; s.Radius > maxRadius+1e-9 {
		t.Errorf("radius %v exceeds half the box diagonal %v", s.Radius, maxRadius)
	}
}
