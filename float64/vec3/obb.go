package vec3

import (
	"math"
)

// obbDegenerateEpsilon is the eigenvalue threshold, relative to the largest
// eigenvalue, below which principal axes are considered degenerate
// (coincident, collinear, coplanar or near-equal eigenvalues).
const obbDegenerateEpsilon = 1e-6

// OBB is an oriented bounding box: a Center, three unit orthonormal Axis
// vectors and the HalfSize extents along each axis.
type OBB struct {
	Center   T
	Axis     [3]T
	HalfSize T
}

// FitOBB fits an oriented bounding box to the given point set along its
// principal directions (PCA of the point covariance matrix).
//
// The returned box contains every input point: HalfSize is enlarged by
// tolerance on each axis, so ContainsPoint reports true for all input points
// even with zero query tolerance. The fit is invariant under translation and
// rotation of the input and under shuffling or duplication of points (up to
// floating point summation order).
//
// Degenerate inputs have a deterministic orientation strategy: coincident
// points yield the canonical axes; collinear points keep the major axis and
// complete the frame with a deterministic perpendicular; coplanar points keep
// the plane normal; near-equal eigenvalues rebuild the ambiguous axes from
// the well-defined one. Every axis is sign-canonicalized so that its
// largest-magnitude component is positive.
func FitOBB(points []T, tolerance float64) OBB {
	obb := OBB{Axis: [3]T{UnitX, UnitY, UnitZ}}
	if len(points) == 0 {
		return obb
	}

	var mean T
	for i := range points {
		mean.Add(&points[i])
	}
	mean.Scale(1 / float64(len(points)))

	var c00, c01, c02, c11, c12, c22 float64
	for i := range points {
		d := Sub(&points[i], &mean)
		c00 += d[0] * d[0]
		c01 += d[0] * d[1]
		c02 += d[0] * d[2]
		c11 += d[1] * d[1]
		c12 += d[1] * d[2]
		c22 += d[2] * d[2]
	}
	cov := [3][3]float64{
		{c00, c01, c02},
		{c01, c11, c12},
		{c02, c12, c22},
	}
	values, vectors := jacobiEigen3(&cov)
	obb.Axis = obbAxes(&values, &vectors)

	// Project all points onto the axes to find the exact extents.
	min := [3]float64{math.MaxFloat64, math.MaxFloat64, math.MaxFloat64}
	max := [3]float64{-math.MaxFloat64, -math.MaxFloat64, -math.MaxFloat64}
	for i := range points {
		d := Sub(&points[i], &mean)
		for a := 0; a < 3; a++ {
			p := Dot(&d, &obb.Axis[a])
			if p < min[a] {
				min[a] = p
			}
			if p > max[a] {
				max[a] = p
			}
		}
	}

	obb.Center = mean
	for a := 0; a < 3; a++ {
		offset := obb.Axis[a].Scaled(0.5 * (min[a] + max[a]))
		obb.Center.Add(&offset)
		obb.HalfSize[a] = 0.5*(max[a]-min[a]) + tolerance
	}
	return obb
}

// ContainsPoint returns whether p lies inside the box, allowing the given
// numerical tolerance on each axis.
func (obb *OBB) ContainsPoint(p *T, tolerance float64) bool {
	d := Sub(p, &obb.Center)
	for a := 0; a < 3; a++ {
		if math.Abs(Dot(&d, &obb.Axis[a])) > obb.HalfSize[a]+tolerance {
			return false
		}
	}
	return true
}

// Corners reconstructs the 8 corners of the box. Corner i uses a negative
// offset along axis a if bit a of i is set, and a positive offset otherwise.
func (obb *OBB) Corners() [8]T {
	var corners [8]T
	for i := 0; i < 8; i++ {
		c := obb.Center
		for a := 0; a < 3; a++ {
			s := float64(1)
			if i&(1<<uint(a)) != 0 {
				s = -1
			}
			offset := obb.Axis[a].Scaled(s * obb.HalfSize[a])
			c.Add(&offset)
		}
		corners[i] = c
	}
	return corners
}

// obbAxes turns raw eigenpairs into a deterministic orthonormal frame,
// applying the degenerate-orientation strategy and sign canonicalization.
func obbAxes(values *[3]float64, vectors *[3]T) [3]T {
	// Sort eigenpairs by eigenvalue, descending (stable for determinism).
	order := [3]int{0, 1, 2}
	for i := 1; i < 3; i++ {
		for j := i; j > 0 && values[order[j-1]] < values[order[j]]; j-- {
			order[j-1], order[j] = order[j], order[j-1]
		}
	}
	var v [3]float64
	var a [3]T
	for i := 0; i < 3; i++ {
		v[i] = values[order[i]]
		a[i] = vectors[order[i]]
	}

	eps := obbDegenerateEpsilon * v[0]
	switch {
	case v[0] <= 0:
		// All points coincide.
		a = [3]T{UnitX, UnitY, UnitZ}
	case v[1] <= eps:
		// Collinear: only the major axis is well defined.
		a[0].Normalize()
		a[1] = obbPerpendicular(&a[0])
		a[2] = Cross(&a[0], &a[1])
	case v[2] <= eps:
		// Coplanar: the plane normal is well defined.
		a[0].Normalize()
		a[1].Normalize()
		a[2] = Cross(&a[0], &a[1])
		a[2].Normalize()
		if v[0]-v[1] <= eps {
			// Near-equal major axes: rebuild both from the normal.
			a[1] = obbPerpendicular(&a[2])
			a[0] = Cross(&a[1], &a[2])
		}
	default:
		a[0].Normalize()
		a[1].Normalize()
		a[2].Normalize()
		if v[0]-v[1] <= eps {
			// Near-equal major axes: rebuild both from the minor axis.
			a[1] = obbPerpendicular(&a[2])
			a[0] = Cross(&a[1], &a[2])
		} else if v[1]-v[2] <= eps {
			// Near-equal minor axes: rebuild both from the major axis.
			a[1] = obbPerpendicular(&a[0])
			a[2] = Cross(&a[0], &a[1])
		}
	}
	for i := range a {
		a[i] = obbCanonicalSign(&a[i])
	}
	return a
}

// obbPerpendicular returns a deterministic unit vector perpendicular to v,
// crossing v with the canonical axis least aligned with it.
func obbPerpendicular(v *T) T {
	ref := UnitX
	if math.Abs(v[1]) < math.Abs(v[0]) {
		ref = UnitY
		if math.Abs(v[2]) < math.Abs(v[1]) {
			ref = UnitZ
		}
	} else if math.Abs(v[2]) < math.Abs(v[0]) {
		ref = UnitZ
	}
	p := Cross(v, &ref)
	return p.Normalized()
}

// obbCanonicalSign flips v so that its largest-magnitude component (lowest
// index on ties) is positive, making the axis orientation deterministic.
func obbCanonicalSign(v *T) T {
	i := 0
	if math.Abs(v[1]) > math.Abs(v[i]) {
		i = 1
	}
	if math.Abs(v[2]) > math.Abs(v[i]) {
		i = 2
	}
	if v[i] < 0 {
		return v.Inverted()
	}
	return *v
}

// jacobiEigen3 computes the eigenvalues and eigenvectors of a symmetric 3x3
// matrix using cyclic Jacobi rotations. vectors[i] is the unit eigenvector
// for values[i].
func jacobiEigen3(m *[3][3]float64) (values [3]float64, vectors [3]T) {
	a := *m
	// v[i][k] is component i of eigenvector k; rotations mix eigenvectors
	// p and q component-wise.
	var v [3][3]float64
	v[0][0], v[1][1], v[2][2] = 1, 1, 1
	for sweep := 0; sweep < 64; sweep++ {
		off := math.Abs(a[0][1]) + math.Abs(a[0][2]) + math.Abs(a[1][2])
		if off == 0 {
			break
		}
		for _, pq := range [][2]int{{0, 1}, {0, 2}, {1, 2}} {
			p, q := pq[0], pq[1]
			apq := a[p][q]
			if apq == 0 {
				continue
			}
			app, aqq := a[p][p], a[q][q]
			theta := (aqq - app) / (2 * apq)
			t := 1 / (math.Abs(theta) + math.Sqrt(theta*theta+1))
			if theta < 0 {
				t = -t
			}
			c := 1 / math.Sqrt(t*t+1)
			s := t * c
			a[p][p] = app - t*apq
			a[q][q] = aqq + t*apq
			a[p][q] = 0
			a[q][p] = 0
			for r := 0; r < 3; r++ {
				if r == p || r == q {
					continue
				}
				arp, arq := a[r][p], a[r][q]
				a[r][p] = c*arp - s*arq
				a[p][r] = a[r][p]
				a[r][q] = s*arp + c*arq
				a[q][r] = a[r][q]
			}
			for i := 0; i < 3; i++ {
				vip, viq := v[i][p], v[i][q]
				v[i][p] = c*vip - s*viq
				v[i][q] = s*vip + c*viq
			}
		}
	}
	values = [3]float64{a[0][0], a[1][1], a[2][2]}
	for k := 0; k < 3; k++ {
		vectors[k] = T{v[0][k], v[1][k], v[2][k]}
	}
	return values, vectors
}
