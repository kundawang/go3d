package vec3

import (
	"errors"
	"fmt"
	stdmath "math"

	math "github.com/ungerik/go3d/fmath"
)

// Sphere is a 3D sphere defined by a center point and a radius.
type Sphere struct {
	Center T
	Radius float32
}

// ErrEmptyPointSet is returned by MinimumBoundingSphere for an empty point set.
var ErrEmptyPointSet = errors.New("vec3: minimum bounding sphere of empty point set")

const (
	// sphereRelativeContainmentTolerance is the relative tolerance used to
	// decide whether a point is contained in a candidate sphere.
	sphereRelativeContainmentTolerance = 1e-4
	// sphereAbsoluteContainmentTolerance is the absolute tolerance used to
	// decide whether a point is contained in a candidate sphere. It avoids
	// a zero tolerance for spheres with a zero radius.
	sphereAbsoluteContainmentTolerance = 1e-6
	// sphereDegenerateTolerance is the relative tolerance used to detect
	// degenerate (collinear or coplanar) point configurations.
	sphereDegenerateTolerance = 1e-6
)

// String formats Sphere as string.
func (s *Sphere) String() string {
	return fmt.Sprint(s.Center.String(), " ", s.Radius)
}

// ContainsPoint returns if a point is contained within the sphere.
func (s *Sphere) ContainsPoint(p *T) bool {
	return SquareDistance(&s.Center, p) <= s.Radius*s.Radius
}

// MinimumBoundingSphere returns the minimum bounding sphere of a finite set
// of 3D points, i.e. the sphere with the smallest possible radius that
// contains all points.
//
// The computation is deterministic: it uses no randomization and does not
// modify the input slice, so repeated calls, and calls with re-ordered or
// duplicated points, return the same result up to floating point precision.
// The returned sphere contains all input points exactly (see ContainsPoint)
// and its radius is minimal up to floating point precision.
//
// An empty point set returns ErrEmptyPointSet. A single point returns a
// sphere with the point as center and a zero radius. Collinear, coplanar
// and near-degenerate point sets are handled with stable fallbacks.
func MinimumBoundingSphere(points []T) (Sphere, error) {
	if len(points) == 0 {
		return Sphere{}, ErrEmptyPointSet
	}
	// Welzl's algorithm, deterministic because the points are processed
	// in input order instead of in a random permutation.
	s := welzl(points, len(points), make([]T, 0, 4))
	// Compensate floating point errors accumulated during the recursion
	// by enlarging the radius so that all points are contained exactly.
	maxDistanceSq := float32(0)
	for i := range points {
		if d := SquareDistance(&s.Center, &points[i]); d > maxDistanceSq {
			maxDistanceSq = d
		}
	}
	s.Radius = math.Sqrt(maxDistanceSq)
	for s.Radius*s.Radius < maxDistanceSq {
		s.Radius = stdmath.Nextafter32(s.Radius, float32(stdmath.Inf(1)))
	}
	return s, nil
}

// welzl implements Welzl's recursive minimum bounding sphere algorithm.
// p is not modified, only the first n points of p are considered and
// boundary holds up to 4 points that are known to lie on the sphere.
func welzl(p []T, n int, boundary []T) Sphere {
	if n == 0 || len(boundary) == 4 {
		return sphereFromBoundary(boundary)
	}
	point := p[n-1]
	s := welzl(p, n-1, boundary)
	if sphereContainsPoint(&s, &point) {
		return s
	}
	return welzl(p, n-1, append(boundary, point))
}

// sphereFromBoundary returns the minimum bounding sphere of up to 4
// boundary points.
func sphereFromBoundary(boundary []T) Sphere {
	switch len(boundary) {
	case 1:
		return Sphere{Center: boundary[0]}
	case 2:
		return sphereFromPair(&boundary[0], &boundary[1])
	case 3:
		return sphereFromTriangle(&boundary[0], &boundary[1], &boundary[2])
	case 4:
		return sphereFromTetrahedron(&boundary[0], &boundary[1], &boundary[2], &boundary[3])
	}
	return Sphere{}
}

// sphereFromPair returns the sphere that has a and b as antipodal points.
func sphereFromPair(a, b *T) Sphere {
	center := Add(a, b)
	center.Scale(0.5)
	return Sphere{Center: center, Radius: Distance(a, b) * 0.5}
}

// sphereFromTriangle returns the sphere through the circumcircle of the
// triangle a, b, c. Collinear (or near collinear) triangles fall back to
// the smallest pair sphere containing all three points.
func sphereFromTriangle(a, b, c *T) Sphere {
	ab := Sub(b, a)
	ac := Sub(c, a)
	n := Cross(&ab, &ac)
	abLenSq := ab.LengthSqr()
	acLenSq := ac.LengthSqr()
	nLenSq := n.LengthSqr()
	if nLenSq <= sphereDegenerateTolerance*sphereDegenerateTolerance*abLenSq*acLenSq {
		return minimumPairSphere([]T{*a, *b, *c})
	}
	// Circumcenter of the triangle:
	// a + (|ab|^2 * (ac x n) + |ac|^2 * (n x ab)) / (2 * |n|^2)
	t1 := Cross(&ac, &n)
	t1.Scale(abLenSq)
	t2 := Cross(&n, &ab)
	t2.Scale(acLenSq)
	t1.Add(&t2)
	t1.Scale(1 / (2 * nLenSq))
	center := Add(a, &t1)
	return Sphere{Center: center, Radius: Distance(&center, a)}
}

// sphereFromTetrahedron returns the circumsphere of the tetrahedron
// a, b, c, d. Coplanar (or near coplanar) point sets fall back to the
// smallest sphere defined by a subset of the points.
func sphereFromTetrahedron(a, b, c, d *T) Sphere {
	r1 := Sub(b, a)
	r2 := Sub(c, a)
	r3 := Sub(d, a)
	c23 := Cross(&r2, &r3)
	det := Dot(&r1, &c23)
	scale := r1.Length() * r2.Length() * r3.Length()
	if math.Abs(det) <= sphereDegenerateTolerance*scale {
		return minimumSphereFromSubsets([]T{*a, *b, *c, *d})
	}
	// Solve the linear system for the circumcenter with Cramer's rule.
	rhs1 := 0.5 * (Dot(b, b) - Dot(a, a))
	rhs2 := 0.5 * (Dot(c, c) - Dot(a, a))
	rhs3 := 0.5 * (Dot(d, d) - Dot(a, a))
	n1 := c23.Scaled(rhs1)
	n2 := Cross(&r3, &r1)
	n2.Scale(rhs2)
	n3 := Cross(&r1, &r2)
	n3.Scale(rhs3)
	n1.Add(&n2).Add(&n3).Scale(1 / det)
	return Sphere{Center: n1, Radius: Distance(&n1, a)}
}

// minimumSphereFromSubsets returns the smallest sphere defined by a pair or
// triple of the given points that contains all of them. It is used as a
// stable fallback for degenerate (coplanar) boundary point sets.
func minimumSphereFromSubsets(points []T) Sphere {
	best := minimumPairSphere(points)
	for i := 0; i < len(points); i++ {
		for j := i + 1; j < len(points); j++ {
			for k := j + 1; k < len(points); k++ {
				s := sphereFromTriangle(&points[i], &points[j], &points[k])
				if s.Radius < best.Radius && sphereContainsAllPoints(&s, points) {
					best = s
				}
			}
		}
	}
	return best
}

// minimumPairSphere returns the smallest sphere around a pair of the given
// points that contains all of them.
func minimumPairSphere(points []T) Sphere {
	best := Sphere{Radius: math.MaxFloat32}
	for i := 0; i < len(points); i++ {
		for j := i + 1; j < len(points); j++ {
			s := sphereFromPair(&points[i], &points[j])
			if s.Radius < best.Radius && sphereContainsAllPoints(&s, points) {
				best = s
			}
		}
	}
	if best.Radius == math.MaxFloat32 && len(points) > 0 {
		// No pair sphere contained all points, which can only happen
		// through floating point error. Fall back to the sphere around
		// the most distant pair of points.
		best = diameterSphere(points)
	}
	return best
}

// diameterSphere returns the sphere around the most distant pair of points.
func diameterSphere(points []T) Sphere {
	a, b := &points[0], &points[0]
	maxDistanceSq := float32(0)
	for i := 0; i < len(points); i++ {
		for j := i + 1; j < len(points); j++ {
			if d := SquareDistance(&points[i], &points[j]); d > maxDistanceSq {
				maxDistanceSq = d
				a, b = &points[i], &points[j]
			}
		}
	}
	return sphereFromPair(a, b)
}

// sphereContainsAllPoints returns if all points are contained in the sphere
// within a small tolerance.
func sphereContainsAllPoints(s *Sphere, points []T) bool {
	for i := range points {
		if !sphereContainsPoint(s, &points[i]) {
			return false
		}
	}
	return true
}

// sphereContainsPoint returns if a point is contained in the sphere within
// a small tolerance.
func sphereContainsPoint(s *Sphere, p *T) bool {
	tolerance := sphereRelativeContainmentTolerance*s.Radius + sphereAbsoluteContainmentTolerance
	d := s.Radius + tolerance
	return SquareDistance(&s.Center, p) <= d*d
}
