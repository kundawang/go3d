package vec3

import (
	"math"
)

// CapsuleHit describes the result of a static capsule vs. axis-aligned box
// collision query.
type CapsuleHit struct {
	// Contact reports whether the capsule touches or overlaps the box.
	// Tangent contact (distance == radius) counts as contact.
	Contact bool
	// Depth is the penetration depth, zero for tangent or non-contacting
	// configurations. Pushing the capsule out of the box by Depth along
	// Normal resolves the overlap.
	Depth float64
	// Normal is the unit contact normal pointing from the box toward the
	// capsule, i.e. the direction the capsule has to be moved to separate
	// both shapes. When the geometry makes the normal non-unique, the
	// smallest axis index wins, so the result is deterministic.
	Normal T
	// BoxPoint is the closest point on the surface of the box.
	BoxPoint T
	// CapsulePoint is the closest point on the surface of the capsule.
	CapsulePoint T
}

// IntersectCapsule performs a static collision query between a capsule
// defined by the endpoints a and b of its axis segment and a radius, and
// the axis-aligned box. A zero-length segment (a == b) degenerates the
// capsule to a sphere.
//
// The query never panics: an invalid box (a Min component greater than the
// corresponding Max component), a negative radius or NaN/Inf input values
// yield a zero CapsuleHit with Contact == false.
//
// Box.Intersects is not modified by this method.
func (box *Box) IntersectCapsule(a, b *T, radius float64) CapsuleHit {
	if radius < 0 ||
		!finiteT(a) || !finiteT(b) || !finiteF(radius) ||
		!finiteT(&box.Min) || !finiteT(&box.Max) {
		return CapsuleHit{}
	}
	for i := 0; i < 3; i++ {
		if box.Min[i] > box.Max[i] {
			return CapsuleHit{}
		}
	}

	// Decide whether the axis segment contains a point strictly inside the
	// box volume. Intersect the open slab intervals (Min,Max) along each
	// axis with the segment parameter [0,1]. A segment merely grazing a
	// face/edge never yields a strict interior point and is handled by the
	// exterior path, which produces a well-defined closest feature.
	tmin, tmax := 0.0, 1.0
	interior := true
	for i := 0; i < 3; i++ {
		d := b[i] - a[i]
		if d == 0 {
			if !(a[i] > box.Min[i] && a[i] < box.Max[i]) {
				interior = false
				break
			}
			continue
		}
		lo := (box.Min[i] - a[i]) / d
		hi := (box.Max[i] - a[i]) / d
		if lo > hi {
			lo, hi = hi, lo
		}
		if lo > tmin {
			tmin = lo
		}
		if hi < tmax {
			tmax = hi
		}
		if tmin >= tmax {
			interior = false
			break
		}
	}
	if interior && tmin < tmax {
		return capsuleBoxInterior(box, a, b, radius)
	}
	return capsuleBoxExterior(box, a, b, radius)
}

// capsuleBoxExterior handles the general case where the capsule axis does
// not pass through the box interior: the closest pair of points between the
// segment and the box is found, and the capsule radius is measured against
// their distance.
func capsuleBoxExterior(box *Box, a, b *T, radius float64) CapsuleHit {
	var hit CapsuleHit

	// distSqr(t) is the squared distance of segment point a+t*(b-a) to the
	// box. Clamping makes it a convex, piecewise quadratic function, so a
	// ternary search converges to the global minimum. Slab crossings are
	// evaluated explicitly as belt-and-braces against flat pieces.
	distSqr := func(t float64) float64 {
		p := Interpolate(a, b, t)
		q := p.Clamped(&box.Min, &box.Max)
		return SquareDistance(&p, &q)
	}

	candidates := make([]float64, 0, 8)
	candidates = append(candidates, 0, 1)
	for i := 0; i < 3; i++ {
		d := b[i] - a[i]
		if d == 0 {
			continue
		}
		for _, plane := range [2]float64{box.Min[i], box.Max[i]} {
			t := (plane - a[i]) / d
			if t > 0 && t < 1 {
				candidates = append(candidates, t)
			}
		}
	}

	// Ternary search on [0,1]; the function is convex so the interval
	// shrinks around a (possibly non-unique) minimum.
	lo, hi := 0.0, 1.0
	for iter := 0; iter < 80; iter++ {
		m1 := lo + (hi-lo)/3
		m2 := hi - (hi-lo)/3
		if distSqr(m1) < distSqr(m2) {
			hi = m2
		} else {
			lo = m1
		}
	}
	candidates = append(candidates, (lo+hi)/2)

	bestSqr := math.Inf(1)
	bestT := 0.0
	for _, t := range candidates {
		sqr := distSqr(t)
		// Ties (e.g. a segment running parallel to a face) resolve toward
		// the segment midpoint so the contact point is deterministic.
		if sqr < bestSqr || (sqr == bestSqr && math.Abs(t-0.5) < math.Abs(bestT-0.5)) {
			bestSqr = sqr
			bestT = t
		}
	}
	// Polishing step: on the final tiny interval the analytic quadratic
	// vertex (same clamp pattern) removes the ternary residue.
	{
		tm := bestT
		span := 1e-7
		t1, t2 := tm-span, tm+span
		if t1 < 0 {
			t1 = 0
		}
		if t2 > 1 {
			t2 = 1
		}
		tmm := (t1 + t2) / 2
		f1, fm, f2 := distSqr(t1), distSqr(tmm), distSqr(t2)
		h := (t2 - t1) / 2
		c2 := (f1 - 2*fm + f2) / (2 * h * h)
		if c2 > 0 {
			c1 := (fm-f1)/h - c2*h
			vertex := tmm - c1/(2*c2)
			if vertex >= 0 && vertex <= 1 {
				sqr := distSqr(vertex)
				if sqr <= bestSqr {
					bestSqr, bestT = sqr, vertex
				}
			}
		}
	}
	// Flat-minimum resolution: when a contiguous sub-interval of the
	// segment shares the global minimum distance (a segment running
	// parallel to a face or alongside an edge), pick the center of that
	// region so the contact location is deterministic and well placed.
	bestSqr = distSqr(bestT)
	tol := 1e-12 * (1 + bestSqr)
	isMin := func(t float64) bool { return distSqr(t) <= bestSqr+tol }
	flatLo, flatHi := 1.0, 0.0
	const samples = 256
	for k := 0; k <= samples; k++ {
		t := float64(k) / samples
		if isMin(t) {
			if t < flatLo {
				flatLo = t
			}
			if t > flatHi {
				flatHi = t
			}
		}
	}
	if flatLo <= flatHi {
		// Refine both boundaries by bisection: just outside the region is
		// a non-minimum, just inside is a minimum.
		bisect := func(out, in float64) float64 {
			for i := 0; i < 50; i++ {
				m := (out + in) / 2
				if isMin(m) {
					in = m
				} else {
					out = m
				}
			}
			return (out + in) / 2
		}
		flatLo = bisect(0, flatLo)
		flatHi = bisect(1, flatHi)
		bestT = (flatLo + flatHi) / 2
		bestSqr = distSqr(bestT)
	}

	segmentPoint := Interpolate(a, b, bestT)
	boxPoint := segmentPoint.Clamped(&box.Min, &box.Max)
	delta := Sub(&segmentPoint, &boxPoint)
	dist := math.Sqrt(bestSqr)

	// Tolerance scaled with the input magnitudes so tangent contact is
	// robust against floating point rounding.
	scale := radius + dist + box.Min.Length() + box.Max.Length() + a.Length() + b.Length()
	epsilon := 1e-12 * scale
	if epsilon == 0 {
		epsilon = 1e-12
	}
	if dist > radius+epsilon {
		// No contact: still report the geometrically closest points and
		// the separation direction.
		hit.BoxPoint = boxPoint
		if dist > 0 {
			hit.Normal = delta.Scaled(1 / dist)
			// Closest point of the capsule surface toward the box.
			offset := hit.Normal.Scaled(-radius)
			hit.CapsulePoint = Add(&segmentPoint, &offset)
		}
		return hit
	}

	hit.Contact = true
	if dist > 0 {
		hit.Normal = delta.Scaled(1 / dist)
		hit.Depth = radius - dist
		if hit.Depth < 0 {
			hit.Depth = 0
		}
	} else {
		// Segment point lies on the box surface, the capsule rests on it.
		hit.Depth = radius
		hit.Normal = faceNormal(box, &boxPoint)
	}
	hit.BoxPoint = boxPoint
	// Closest point on the capsule surface, exactly Depth away from the box.
	offset := hit.Normal.Scaled(-radius)
	hit.CapsulePoint = Add(&segmentPoint, &offset)
	return hit
}

// capsuleBoxInterior handles the case where the capsule axis passes through
// the box volume, so no unique closest feature exists. The capsule extent on
// every axis is [lo-r, hi+r]; wherever it crosses a box face the positive
// overlap is a candidate minimum translation vector. The smallest one is
// chosen. Ties resolve to the smallest axis index and then the low face, so
// the result is deterministic when the normal is non-unique.
func capsuleBoxInterior(box *Box, a, b *T, radius float64) CapsuleHit {
	var hit CapsuleHit

	bestDepth := math.Inf(1)
	bestAxis := 0
	bestSign := -1.0
	bestAxisPoint := T{}
	for i := 0; i < 3; i++ {
		lo, hi := a[i], b[i]
		if lo > hi {
			lo, hi = hi, lo
		}
		// Overlap through the low face (outward direction -axis i):
		// Min - (lo - r), only positive means a real overlap.
		if d := box.Min[i] - (lo - radius); d > 0 && d < bestDepth {
			bestDepth = d
			bestAxis = i
			bestSign = -1
			bestAxisPoint = minEndpoint(a, b, i)
		}
		// Overlap through the high face (outward direction +axis i):
		// (hi + r) - Max.
		if d := (hi + radius) - box.Max[i]; d > 0 && d < bestDepth {
			bestDepth = d
			bestAxis = i
			bestSign = +1
			bestAxisPoint = maxEndpoint(a, b, i)
		}
	}

	hit.Contact = true
	if bestDepth == math.Inf(1) {
		// No face overlap despite an interior axis: treat as tangent.
		bestDepth = 0
	}
	hit.Depth = bestDepth
	hit.Normal = T{}
	hit.Normal[bestAxis] = bestSign

	// Box surface point: the face plane with the other coordinates clamped
	// to the face rectangle.
	facePlane := box.Min[bestAxis]
	if bestSign > 0 {
		facePlane = box.Max[bestAxis]
	}
	boxPoint := bestAxisPoint.Clamped(&box.Min, &box.Max)
	boxPoint[bestAxis] = facePlane
	hit.BoxPoint = boxPoint
	// Capsule surface point: exactly Depth away from the box point along the
	// separating normal, giving the uniform invariant
	// CapsulePoint - BoxPoint == Depth * Normal for both code paths.
	depthOffset := hit.Normal.Scaled(hit.Depth)
	hit.CapsulePoint = Add(&boxPoint, &depthOffset)
	return hit
}

// minEndpoint returns the endpoint a/b with the smaller component on axis.
func minEndpoint(a, b *T, axis int) T {
	if a[axis] <= b[axis] {
		return *a
	}
	return *b
}

// maxEndpoint returns the endpoint a/b with the larger component on axis.
func maxEndpoint(a, b *T, axis int) T {
	if a[axis] >= b[axis] {
		return *a
	}
	return *b
}

// faceNormal returns a deterministic outward normal for a point located on
// the surface of the box, preferring the smallest axis index when several
// faces meet (edges and corners).
func faceNormal(box *Box, p *T) T {
	for i := 0; i < 3; i++ {
		if p[i] == box.Min[i] {
			n := T{}
			n[i] = -1
			return n
		}
		if p[i] == box.Max[i] {
			n := T{}
			n[i] = +1
			return n
		}
	}
	return T{1, 0, 0}
}

func finiteF(v float64) bool {
	return !math.IsNaN(v) && !math.IsInf(v, 0)
}

func finiteT(v *T) bool {
	for i := 0; i < 3; i++ {
		if !finiteF(v[i]) {
			return false
		}
	}
	return true
}
