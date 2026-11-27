package quaternion_test

import (
	. "github.com/ungerik/go3d/quaternion"
	"math"
	"math/rand"
	"testing"

	math32 "github.com/chewxy/math32"
	quat64 "github.com/ungerik/go3d/float64/quaternion"
	"github.com/ungerik/go3d/vec3"
)

// uid127SameOrientation checks whether two unit quaternions represent the same
// rotation, which is the case when they are equal or exact negatives of
// each other (the double cover of SO(3)).
func uid127SameOrientation(a, b T, tolerance float32) bool {
	equal := true
	negated := true
	for i := 0; i < 4; i++ {
		if math32.Abs(a[i]-b[i]) > tolerance {
			equal = false
		}
		if math32.Abs(a[i]+b[i]) > tolerance {
			negated = false
		}
	}
	return equal || negated
}

// uid127ZAngle returns the rotation angle of a quaternion that is known to rotate
// around the Z axis, computed directly from the scalar component so that
// AxisAngle()'s float32 precision behavior does not affect the assertion.
func uid127ZAngle(q T) float32 {
	return 2 * math32.Atan2(q[2], q[3])
}

func TestUID127WeightedAverageTwoQuats45Degrees(t *testing.T) {
	q0 := FromZAxisAngle(0)
	q90 := FromZAxisAngle(math32.Pi / 2)

	result, err := WeightedAverage([]T{q0, q90}, []float32{1, 1})
	if err != nil {
		t.Fatalf("WeightedAverage() error = %v", err)
	}
	if !result.IsUnitQuat(1e-6) {
		t.Fatalf("result is not a unit quaternion: %v (norm=%v)", result, result.Norm())
	}
	if angle := uid127ZAngle(result); math32.Abs(angle-math32.Pi/4) > 1e-5 {
		t.Errorf("average angle = %v, want %v (45 degrees)", angle, math32.Pi/4)
	}
	if math32.Abs(result[0]) > 1e-6 || math32.Abs(result[1]) > 1e-6 {
		t.Errorf("average should rotate around Z, got %v", result)
	}
}

func TestUID127WeightedAverageNegationInvariance(t *testing.T) {
	q0 := FromZAxisAngle(0)
	q90 := FromZAxisAngle(math32.Pi / 2)

	base, err := WeightedAverage([]T{q0, q90}, []float32{1, 1})
	if err != nil {
		t.Fatalf("WeightedAverage() error = %v", err)
	}
	// Negating any input quaternion leaves the represented orientation,
	// and therefore the average, unchanged.
	for _, neg := range [][]int{{0}, {1}, {0, 1}} {
		quats := []T{q0, q90}
		for _, i := range neg {
			quats[i] = quats[i].Negated()
		}
		result, err := WeightedAverage(quats, []float32{1, 1})
		if err != nil {
			t.Fatalf("WeightedAverage(negated %v) error = %v", neg, err)
		}
		if !uid127SameOrientation(base, result, 1e-6) {
			t.Errorf("negating inputs %v changed orientation: base=%v result=%v", neg, base, result)
		}
		// Canonical sign makes the quaternions bit-for-bit equivalent too.
		if result != base {
			t.Errorf("negating inputs %v changed canonical quaternion: base=%v result=%v", neg, base, result)
		}
	}
}

func TestUID127WeightedAverageSingleEffectivePose(t *testing.T) {
	q := FromAxisAngle(&vec3.T{1, 2, 3}, 0.7)

	// A single quaternion with positive weight is returned as itself.
	result, err := WeightedAverage([]T{q}, []float32{2.5})
	if err != nil {
		t.Fatalf("WeightedAverage() error = %v", err)
	}
	if !uid127SameOrientation(result, q, 1e-6) {
		t.Errorf("single pose average = %v, want %v", result, q)
	}

	// Zero-weight companions are ignored, so the only weighted pose wins.
	other := FromZAxisAngle(math32.Pi)
	result, err = WeightedAverage([]T{other, q}, []float32{0, 1})
	if err != nil {
		t.Fatalf("WeightedAverage() error = %v", err)
	}
	if !uid127SameOrientation(result, q, 1e-6) {
		t.Errorf("zero-weight companion changed result: got %v, want %v", result, q)
	}

	// Same orientation presented with both signs still gives the pose.
	result, err = WeightedAverage([]T{q, q.Negated()}, []float32{1, 1})
	if err != nil {
		t.Fatalf("WeightedAverage() error = %v", err)
	}
	if !uid127SameOrientation(result, q, 1e-6) {
		t.Errorf("q/-q average = %v, want orientation %v", result, q)
	}
}

func TestUID127WeightedAverageZeroWeights(t *testing.T) {
	q0 := FromZAxisAngle(0)
	q90 := FromZAxisAngle(math32.Pi / 2)

	result, err := WeightedAverage([]T{q0, q90}, []float32{0, 0})
	if err != ErrAllWeightsZero {
		t.Errorf("all-zero weights: err = %v, want ErrAllWeightsZero", err)
	}
	if result != Zero {
		t.Errorf("all-zero weights: result = %v, want Zero", result)
	}

	if _, err := WeightedAverage(nil, nil); err != ErrEmptyInput {
		t.Errorf("empty input: err = %v, want ErrEmptyInput", err)
	}
}

func TestUID127WeightedAverageOppositeOrientations(t *testing.T) {
	// Two genuinely opposite orientations (180 degrees apart in SO(3)).
	// The average is not unique physically, but it must still be a
	// deterministic valid unit quaternion instead of failing or producing
	// NaN/zero the way a cancelling vector sum would.
	qx := FromXAxisAngle(math32.Pi)
	qz := FromZAxisAngle(math32.Pi)

	result, err := WeightedAverage([]T{qx, qz}, []float32{1, 1})
	if err != nil {
		t.Fatalf("WeightedAverage() error = %v", err)
	}
	if !result.IsUnitQuat(1e-6) {
		t.Fatalf("opposite orientations gave non-unit result: %v (norm=%v)", result, result.Norm())
	}
	for i := range result {
		if math32.IsNaN(result[i]) || math32.IsInf(result[i], 0) {
			t.Fatalf("opposite orientations gave non-finite result: %v", result)
		}
	}
	again, err := WeightedAverage([]T{qz, qx}, []float32{1, 1})
	if err != nil {
		t.Fatalf("WeightedAverage() error = %v", err)
	}
	if again != result {
		t.Errorf("opposite orientations not deterministic: %v vs %v", result, again)
	}
}

func TestUID127WeightedAverageNear180Degrees(t *testing.T) {
	// Two orientations separated by nearly 180 degrees: a plain weighted
	// vector sum would nearly cancel here.
	small := FromZAxisAngle(0.01)
	large := FromZAxisAngle(math32.Pi - 0.01)

	result, err := WeightedAverage([]T{small, large}, []float32{1, 1})
	if err != nil {
		t.Fatalf("WeightedAverage() error = %v", err)
	}
	if !result.IsUnitQuat(1e-6) {
		t.Fatalf("near-180 spread gave non-unit result: %v (norm=%v)", result, result.Norm())
	}
	// Equal weights must land on the geodesic midpoint around 90 degrees.
	if angle := uid127ZAngle(result); math32.Abs(angle-math32.Pi/2) > 1e-4 {
		t.Errorf("near-180 midpoint angle = %v, want ~%v", angle, math32.Pi/2)
	}
	if math32.Abs(result[0]) > 1e-5 || math32.Abs(result[1]) > 1e-5 {
		t.Errorf("near-180 midpoint should stay around Z axis, got %v", result)
	}
}

func TestUID127WeightedAverageThreeQuatsUnequalWeights(t *testing.T) {
	// 0, 90 and 180 degrees around Z with weights 1, 1, 2.
	// The dominant eigenvector of the resulting 2x2 block represents
	// a 135 degrees rotation around Z.
	q0 := FromZAxisAngle(0)
	q90 := FromZAxisAngle(math32.Pi / 2)
	q180 := FromZAxisAngle(math32.Pi)

	result, err := WeightedAverage([]T{q0, q90, q180}, []float32{1, 1, 2})
	if err != nil {
		t.Fatalf("WeightedAverage() error = %v", err)
	}
	if !result.IsUnitQuat(1e-6) {
		t.Fatalf("result is not a unit quaternion: %v (norm=%v)", result, result.Norm())
	}
	if angle := uid127ZAngle(result); math32.Abs(angle-3*math32.Pi/4) > 1e-5 {
		t.Errorf("three-pose angle = %v, want %v (135 degrees)", angle, 3*math32.Pi/4)
	}
	if math32.Abs(result[0]) > 1e-6 || math32.Abs(result[1]) > 1e-6 {
		t.Errorf("three-pose average should rotate around Z, got %v", result)
	}

	// Larger weight on q0 must pull the two-pose average toward 0 degrees.
	biased, err := WeightedAverage([]T{q0, q90}, []float32{3, 1})
	if err != nil {
		t.Fatalf("WeightedAverage() error = %v", err)
	}
	if a := uid127ZAngle(biased); a <= 0 || a >= math32.Pi/4 {
		t.Errorf("weight 3:1 average angle = %v, want between 0 and 45 degrees", a)
	}
}

func TestUID127WeightedAverageOrderAndSignInvariance(t *testing.T) {
	quats := []T{
		FromAxisAngle(&vec3.UnitX, 0.4),
		FromAxisAngle(&vec3.UnitY, -1.1),
		FromAxisAngle(&vec3.UnitZ, 1.7),
		FromAxisAngle(&vec3.T{1, 1, 0}, 0.9),
	}
	weights := []float32{0.5, 2, 1.25, 0.75}

	base, err := WeightedAverage(quats, weights)
	if err != nil {
		t.Fatalf("WeightedAverage() error = %v", err)
	}

	permutedQuats := []T{quats[3], quats[1], quats[0], quats[2]}
	permutedWeights := []float32{weights[3], weights[1], weights[0], weights[2]}
	permuted, err := WeightedAverage(permutedQuats, permutedWeights)
	if err != nil {
		t.Fatalf("WeightedAverage() error = %v", err)
	}
	// float32 summation order introduces a few ULP of roundoff, so compare
	// the represented orientation rather than individual bits.
	if !uid127SameOrientation(base, permuted, 1e-6) {
		t.Errorf("permutation changed result: base=%v permuted=%v", base, permuted)
	}

	// Negate a subset of quaternions at arbitrary positions.
	signed := make([]T, len(quats))
	copy(signed, quats)
	signed[0] = signed[0].Negated()
	signed[2] = signed[2].Negated()
	result, err := WeightedAverage(signed, weights)
	if err != nil {
		t.Fatalf("WeightedAverage() error = %v", err)
	}
	if result != base {
		t.Errorf("negating inputs changed result: base=%v result=%v", base, result)
	}
}

func TestUID127WeightedAverageNormalizesNonUnitInput(t *testing.T) {
	q := FromAxisAngle(&vec3.UnitZ, 0.6)
	scaled := T{q[0] * 3, q[1] * 3, q[2] * 3, q[3] * 3}

	result, err := WeightedAverage([]T{scaled}, []float32{1})
	if err != nil {
		t.Fatalf("WeightedAverage() error = %v", err)
	}
	if !uid127SameOrientation(result, q, 1e-6) {
		t.Errorf("non-unit input: got %v, want orientation %v", result, q)
	}
}

func TestUID127WeightedAverageExtremeButFiniteScales(t *testing.T) {
	q := FromAxisAngle(&vec3.T{1, 2, 3}, 0.7)
	for _, scale := range []float32{1e-30, 1e30} {
		scaled := T{q[0] * scale, q[1] * scale, q[2] * scale, q[3] * scale}
		result, err := WeightedAverage([]T{scaled}, []float32{1})
		if err != nil {
			t.Fatalf("scale %v: WeightedAverage() error = %v", scale, err)
		}
		if !uid127SameOrientation(result, q, 1e-5) {
			t.Errorf("scale %v changed orientation: got %v, want %v", scale, result, q)
		}
	}

	quats := []T{
		FromZAxisAngle(0),
		FromZAxisAngle(math32.Pi / 3),
		FromZAxisAngle(math32.Pi / 2),
	}
	baseWeights := []float32{1, 2, 3}
	base, err := WeightedAverage(quats, baseWeights)
	if err != nil {
		t.Fatal(err)
	}
	small := float32(math.SmallestNonzeroFloat32)
	huge := float32(math32.MaxFloat32)
	for _, weights := range [][]float32{
		{small, 2 * small, 3 * small},
		{huge / 3, 2 * (huge / 3), huge},
	} {
		result, err := WeightedAverage(quats, weights)
		if err != nil {
			t.Fatalf("weights %v: %v", weights, err)
		}
		if !uid127SameOrientation(base, result, 1e-5) {
			t.Errorf("weight scale changed result: base=%v got=%v weights=%v", base, result, weights)
		}
	}
}

func TestUID127WeightedAverageInvalidInputs(t *testing.T) {
	q := FromZAxisAngle(0.5)

	tests := []struct {
		name    string
		quats   []T
		weights []float32
		wantErr error
	}{
		{"empty", nil, nil, ErrEmptyInput},
		{"length mismatch", []T{q, q}, []float32{1}, ErrInputLengthMismatch},
		{"negative weight", []T{q}, []float32{-1}, ErrNegativeWeight},
		{"NaN weight", []T{q}, []float32{float32(math.NaN())}, ErrNonFiniteInput},
		{"Inf weight", []T{q}, []float32{float32(math.Inf(1))}, ErrNonFiniteInput},
		{"NaN quaternion", []T{{0, 0, 0, float32(math.NaN())}}, []float32{1}, ErrNonFiniteInput},
		{"Inf quaternion", []T{{0, 0, 0, float32(math.Inf(1))}}, []float32{1}, ErrNonFiniteInput},
		{"zero quaternion positive weight", []T{Zero}, []float32{1}, ErrZeroQuaternion},
		{"all weights zero", []T{q, q}, []float32{0, 0}, ErrAllWeightsZero},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result, err := WeightedAverage(tt.quats, tt.weights)
			if err != tt.wantErr {
				t.Errorf("err = %v, want %v", err, tt.wantErr)
			}
			if err != nil && result != Zero {
				t.Errorf("result on error = %v, want Zero", result)
			}
		})
	}

	// Non-finite inputs are validated even when their weight is zero.
	if _, err := WeightedAverage(
		[]T{q, {0, 0, 0, float32(math.NaN())}},
		[]float32{1, 0},
	); err != ErrNonFiniteInput {
		t.Errorf("zero-weight NaN quaternion: err = %v, want ErrNonFiniteInput", err)
	}
}

func TestUID127WeightedAverageFloat32Float64Correspondence(t *testing.T) {
	quats32 := []T{
		FromAxisAngle(&vec3.UnitX, 0.3),
		FromAxisAngle(&vec3.UnitY, 1.2),
		FromAxisAngle(&vec3.UnitZ, -0.8),
	}
	weights32 := []float32{1.5, 0.25, 2}

	quats64 := make([]quat64.T, len(quats32))
	weights64 := make([]float64, len(weights32))
	for i := range quats32 {
		for k := 0; k < 4; k++ {
			quats64[i][k] = float64(quats32[i][k])
		}
		weights64[i] = float64(weights32[i])
	}

	r32, err := WeightedAverage(quats32, weights32)
	if err != nil {
		t.Fatalf("float32 WeightedAverage() error = %v", err)
	}
	r64, err := quat64.WeightedAverage(quats64, weights64)
	if err != nil {
		t.Fatalf("float64 WeightedAverage() error = %v", err)
	}

	const tolerance = 2e-6
	for k := 0; k < 4; k++ {
		if d := math32.Abs(r32[k] - float32(r64[k])); d > tolerance {
			t.Errorf("component %d: float32=%v float64=%v diff=%v", k, r32[k], r64[k], d)
		}
	}
}

func BenchmarkUID127WeightedAverage(b *testing.B) {
	quats := []T{
		FromAxisAngle(&vec3.UnitX, 0.4),
		FromAxisAngle(&vec3.UnitY, -1.1),
		FromAxisAngle(&vec3.UnitZ, 1.7),
		FromAxisAngle(&vec3.T{1, 1, 0}, 0.9),
	}
	weights := []float32{0.5, 2, 1.25, 0.75}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		if _, err := WeightedAverage(quats, weights); err != nil {
			b.Fatal(err)
		}
	}
}

// uid127RandomUnitQuat builds a uniformly-ish random unit quaternion.
func uid127RandomUnitQuat(rng *rand.Rand) T {
	q := T{
		rng.Float32()*2 - 1,
		rng.Float32()*2 - 1,
		rng.Float32()*2 - 1,
		rng.Float32()*2 - 1,
	}
	return q.Normalized()
}

// uid127PowerIterationDominant independently approximates the dominant eigenvector
// of a symmetric 4x4 matrix by repeated matrix-vector multiplication.
func uid127PowerIterationDominant(m *[4][4]float32) [4]float32 {
	v := [4]float32{1, 0.5, -0.25, 0.125}
	for iter := 0; iter < 200; iter++ {
		var next [4]float32
		var norm float32
		for r := 0; r < 4; r++ {
			for c := 0; c < 4; c++ {
				next[r] += m[r][c] * v[c]
			}
			norm += next[r] * next[r]
		}
		if norm == 0 {
			return next
		}
		ool := 1 / math32.Sqrt(norm)
		for r := range next {
			next[r] *= ool
		}
		v = next
	}
	return v
}

func TestUID127WeightedAverageRandomizedStress(t *testing.T) {
	rng := rand.New(rand.NewSource(42))
	for iter := 0; iter < 200; iter++ {
		n := 1 + rng.Intn(6)
		quats := make([]T, n)
		weights := make([]float32, n)
		var m [4][4]float32
		for i := range quats {
			quats[i] = uid127RandomUnitQuat(rng)
			weights[i] = rng.Float32() * 10
			for r := 0; r < 4; r++ {
				for c := r; c < 4; c++ {
					m[r][c] += weights[i] * quats[i][r] * quats[i][c]
				}
			}
		}
		for r := 1; r < 4; r++ {
			for c := 0; c < r; c++ {
				m[r][c] = m[c][r]
			}
		}

		result, err := WeightedAverage(quats, weights)
		if err != nil {
			t.Fatalf("iter %d: unexpected error %v", iter, err)
		}
		ev := [4]float32(result)

		// 1. The result is a unit quaternion.
		var norm float32
		for _, x := range ev {
			norm += x * x
		}
		if math32.Abs(norm-1) > 1e-5 {
			t.Fatalf("iter %d: result norm = %v", iter, norm)
		}

		// 2. It is an eigenvector of the outer-product matrix: M*ev = lambda*ev.
		//    Compute M*ev once per row, then both the Rayleigh quotient and
		//    the eigenvector residual ||M*ev - lambda*ev||.
		var mEv [4]float32
		var rqResult float32
		for r := 0; r < 4; r++ {
			for c := 0; c < 4; c++ {
				mEv[r] += m[r][c] * ev[c]
			}
			rqResult += ev[r] * mEv[r]
		}
		var resid float32
		for r := 0; r < 4; r++ {
			d := mEv[r] - rqResult*ev[r]
			resid += d * d
		}
		if math32.Sqrt(resid) > 1e-3 {
			t.Fatalf("iter %d: eigenvector residual = %v", iter, math32.Sqrt(resid))
		}

		// 3. Dominance: its Rayleigh quotient must match the value found
		//    independently by power iteration (a lower bound on lambda max).
		pi := uid127PowerIterationDominant(&m)
		var rqPower float32
		for r := 0; r < 4; r++ {
			var mPi float32
			for c := 0; c < 4; c++ {
				mPi += m[r][c] * pi[c]
			}
			rqPower += pi[r] * mPi
		}
		if rqResult+1e-3 < rqPower {
			t.Fatalf("iter %d: not dominant, rq=%v power=%v", iter, rqResult, rqPower)
		}

		// 4. Permuting and negating inputs keeps the same orientation.
		order := rng.Perm(n)
		permQuats := make([]T, n)
		permWeights := make([]float32, n)
		for i, j := range order {
			permQuats[i] = quats[j]
			permWeights[i] = weights[j]
			if rng.Intn(2) == 0 {
				permQuats[i] = permQuats[i].Negated()
			}
		}
		again, err := WeightedAverage(permQuats, permWeights)
		if err != nil {
			t.Fatalf("iter %d: permuted error %v", iter, err)
		}
		if !uid127SameOrientation(result, again, 1e-5) {
			t.Fatalf("iter %d: permutation/negation changed orientation: %v vs %v", iter, result, again)
		}
	}
}
