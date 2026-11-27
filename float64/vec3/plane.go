package vec3

// Plane represents an oriented plane in 3D.
//
// For a point p the signed distance is Dot(Normal, p) + Offset.
// Points with a non-negative signed distance lie on the positive side
// of the plane, which is the side that ClipPolygon keeps.
//
// The normal does not have to be unit length. ClipPolygon normalizes it
// internally, so the scale of Normal and Offset does not influence the
// result as long as both are scaled by the same factor.
type Plane struct {
	Normal T
	Offset float64
}

// NewPlane returns the plane through point with the given normal.
//
// With a non-zero normal, points p with Dot(normal, p-point) >= 0 are on
// the positive (kept) side. The normal is stored as given and is not
// normalized.
func NewPlane(point, normal *T) Plane {
	return Plane{
		Normal: *normal,
		Offset: -Dot(normal, point),
	}
}

// SignedDistance returns Dot(plane.Normal, point) + plane.Offset.
//
// For a unit length normal this is the Euclidean distance to the plane.
func (plane *Plane) SignedDistance(point *T) float64 {
	return Dot(&plane.Normal, point) + plane.Offset
}

// ContainsPoint reports whether point lies on the positive side of the
// plane or exactly on the plane. It does not use any tolerance, so a point
// that is only numerically close to the plane may be reported as outside.
func (plane *Plane) ContainsPoint(point *T) bool {
	return plane.SignedDistance(point) >= 0
}
