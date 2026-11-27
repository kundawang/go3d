package quaternion

import (
	"errors"

	math "github.com/chewxy/math32"
)

// Errors returned by WeightedAverage.
var (
	// ErrEmptyInput is returned when no quaternions are provided.
	ErrEmptyInput = errors.New("quaternion.WeightedAverage: no quaternions provided")

	// ErrInputLengthMismatch is returned when the number of quaternions
	// does not match the number of weights.
	ErrInputLengthMismatch = errors.New("quaternion.WeightedAverage: quaternion and weight counts differ")

	// ErrNegativeWeight is returned when at least one weight is negative.
	ErrNegativeWeight = errors.New("quaternion.WeightedAverage: weights must be non-negative")

	// ErrNonFiniteInput is returned when a quaternion component or a weight
	// is NaN or Inf. Entries with zero weight are validated as well.
	// It is also returned when finite inputs make the weighted sum overflow.
	ErrNonFiniteInput = errors.New("quaternion.WeightedAverage: quaternion or weight contains NaN or Inf, or the weighted sum overflows")

	// ErrZeroQuaternion is returned when a quaternion with positive weight
	// has zero norm and therefore represents no orientation.
	ErrZeroQuaternion = errors.New("quaternion.WeightedAverage: quaternion with positive weight has zero norm")

	// ErrAllWeightsZero is returned when every weight is zero, so no
	// average orientation can be derived.
	ErrAllWeightsZero = errors.New("quaternion.WeightedAverage: total weight is zero")
)

// WeightedAverage returns the weighted average orientation of the given
// unit quaternions as a unit quaternion.
//
// It implements the averaging method of Markley et al. (2007):
// the orientation is the normalized eigenvector belonging to the largest
// eigenvalue of M = sum(weight[i] * quat[i] * quat[i]^T).
//
// Properties:
//   - The result does not depend on the order of the input pairs.
//   - Quaternions q and -q represent the same orientation, so negating any
//     input quaternion does not change the result.
//   - Inputs spread by nearly 180 degrees are handled without the
//     cancellation that a plain weighted vector sum would suffer.
//
// Quaternions with a zero weight are ignored, but all inputs are still
// validated. Non-unit quaternions with positive weight are normalized,
// while zero-norm quaternions cause an ErrZeroQuaternion error.
// NaN/Inf components or weights, negative weights, empty input, mismatched
// slice lengths and an all-zero weight sum return a Zero quaternion together
// with the corresponding error.
func WeightedAverage(quats []T, weights []float32) (T, error) {
	n := len(quats)
	if n == 0 {
		return Zero, ErrEmptyInput
	}
	if len(weights) != n {
		return Zero, ErrInputLengthMismatch
	}

	// Validate first and find the largest positive weight. Dividing all
	// weights by that maximum below preserves the eigenproblem while making
	// the result independent of the weights' common scale and avoiding
	// overflow/underflow for very large or small finite weights.
	var maxWeight float32
	for i := 0; i < n; i++ {
		w := weights[i]
		if math.IsNaN(w) || math.IsInf(w, 0) {
			return Zero, ErrNonFiniteInput
		}
		if w < 0 {
			return Zero, ErrNegativeWeight
		}
		for k := 0; k < 4; k++ {
			if math.IsNaN(quats[i][k]) || math.IsInf(quats[i][k], 0) {
				return Zero, ErrNonFiniteInput
			}
		}
		if w > maxWeight {
			maxWeight = w
		}
	}
	if maxWeight == 0 {
		return Zero, ErrAllWeightsZero
	}

	var m [4][4]float32
	for i := 0; i < n; i++ {
		w := weights[i]
		if w == 0 {
			continue
		}
		q, ok := normalizedAverageQuaternion(quats[i])
		if !ok {
			return Zero, ErrZeroQuaternion
		}
		normalizedWeight := w / maxWeight
		for r := 0; r < 4; r++ {
			for c := r; c < 4; c++ {
				m[r][c] += normalizedWeight * q[r] * q[c]
			}
		}
	}
	// Reject overflow from an extreme number of terms so that no NaN/Inf
	// can reach the Jacobi iteration from otherwise finite inputs.
	for r := 0; r < 4; r++ {
		for c := r; c < 4; c++ {
			if math.IsNaN(m[r][c]) || math.IsInf(m[r][c], 0) {
				return Zero, ErrNonFiniteInput
			}
		}
	}
	// Mirror the upper triangle into the lower one for Jacobi iteration.
	for r := 1; r < 4; r++ {
		for c := 0; c < r; c++ {
			m[r][c] = m[c][r]
		}
	}

	ev := dominantEigenvector(&m)

	// Every orientation has two quaternion representations (ev and -ev).
	// Fix a deterministic sign: the first non-zero component is positive.
	for r := 0; r < 4; r++ {
		if ev[r] != 0 {
			if ev[r] < 0 {
				for k := range ev {
					ev[k] = -ev[k]
				}
			}
			break
		}
	}

	return T(ev), nil
}

// normalizedAverageQuaternion normalizes q for averaging. Scaling by the
// largest absolute component first keeps extremely small finite quaternions
// stable when Norm() would underflow or fall below the package Epsilon.
func normalizedAverageQuaternion(q T) (T, bool) {
	maxAbs := math.Abs(q[0])
	for i := 1; i < len(q); i++ {
		if a := math.Abs(q[i]); a > maxAbs {
			maxAbs = a
		}
	}
	if maxAbs == 0 {
		return q, false
	}

	var scaled T
	var normSquared float32
	for i := range q {
		scaled[i] = q[i] / maxAbs
		normSquared += scaled[i] * scaled[i]
	}
	inverseNorm := 1 / math.Sqrt(normSquared)
	for i := range scaled {
		scaled[i] *= inverseNorm
	}
	return scaled, true
}

// dominantEigenvector returns a unit eigenvector belonging to the largest
// eigenvalue of the symmetric 4x4 matrix a, computed by cyclic Jacobi
// eigenvalue iteration. The fixed pivot order makes the result deterministic
// even when eigenvalues are repeated. The matrix is overwritten during
// iteration, with its diagonal holding the eigenvalues afterwards.
func dominantEigenvector(a *[4][4]float32) [4]float32 {
	// v accumulates the Jacobi rotations and holds the eigenvectors as
	// columns when the iteration converges.
	var v [4][4]float32
	for i := 0; i < 4; i++ {
		v[i][i] = 1
	}

	const maxSweeps = 30
	for sweep := 0; sweep < maxSweeps; sweep++ {
		converged := true
		for p := 0; p < 3; p++ {
			for q := p + 1; q < 4; q++ {
				apq := a[p][q]
				if math.Abs(apq) < 1e-7*(1+math.Abs(a[p][p])+math.Abs(a[q][q])) {
					continue
				}
				converged = false

				// tan(theta) solving the zeroing equation for a[p][q],
				// stabilized so that |t| <= 1.
				tau := (a[q][q] - a[p][p]) / (2 * apq)
				var t float32
				if tau >= 0 {
					t = 1 / (tau + math.Sqrt(1+tau*tau))
				} else {
					t = -1 / (-tau + math.Sqrt(1+tau*tau))
				}
				c := 1 / math.Sqrt(1+t*t)
				s := t * c

				app := a[p][p]
				aqq := a[q][q]
				a[p][p] = c*c*app - 2*s*c*apq + s*s*aqq
				a[q][q] = s*s*app + 2*s*c*apq + c*c*aqq
				a[p][q] = 0
				a[q][p] = 0
				for r := 0; r < 4; r++ {
					if r == p || r == q {
						continue
					}
					arp := a[r][p]
					arq := a[r][q]
					a[r][p] = c*arp - s*arq
					a[p][r] = a[r][p]
					a[r][q] = s*arp + c*arq
					a[q][r] = a[r][q]
				}
				for r := 0; r < 4; r++ {
					vrp := v[r][p]
					vrq := v[r][q]
					v[r][p] = c*vrp - s*vrq
					v[r][q] = s*vrp + c*vrq
				}
			}
		}
		if converged {
			break
		}
	}

	// The largest diagonal entry is the largest eigenvalue; its column in v
	// is the corresponding eigenvector. Strict comparison keeps the choice
	// deterministic for repeated eigenvalues.
	best := 0
	for i := 1; i < 4; i++ {
		if a[i][i] > a[best][best] {
			best = i
		}
	}

	var ev [4]float32
	var norm float32
	for i := 0; i < 4; i++ {
		ev[i] = v[i][best]
		norm += ev[i] * ev[i]
	}
	ool := 1 / math.Sqrt(norm)
	for i := range ev {
		ev[i] *= ool
	}
	return ev
}
