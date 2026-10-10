"""Optimal two-qubit synthesis (KAK / Cartan decomposition).

Every two-qubit unitary can be written as

    U = e^{iγ} (A1 ⊗ A0) · exp(i(a XX + b YY + c ZZ)) · (B1 ⊗ B0)

with single-qubit ``A``s and ``B``s and a point ``(a, b, c)`` of the Weyl
chamber ``π/4 ≥ a ≥ b ≥ |c|`` (with ``c ≥ 0`` when ``a = π/4``). Two gates
are equal up to single-qubit gates exactly when their points coincide, and the
point fixes how many uses of a native two-qubit gate a decomposition needs:

* CX-like gates (cx, cz, cy, ch, ecr, ...): 0 for (0, 0, 0), 1 for (π/4, 0, 0),
  2 when c = 0 and 3 otherwise. These counts are optimal.
* Ising rotations (rxx, ryy, rzz, rzx) with a free angle: one per non-zero
  coordinate, at most three.

The decomposition writes ``U`` as the canonical gate of a fixed template
circuit with the right point, sandwiched between single-qubit gates computed
from the KAK decompositions of ``U`` and of the template.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Collection, Sequence
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from qevedo.compiler.synthesis.gates import SWAP, H, S, X, Y, Z, gate_matrix, rx, rz
from qevedo.compiler.synthesis.gates import SX as _SX

__all__ = ["KAK", "TwoQubitBasis", "best_2q_bases", "kak", "synthesize_2q", "two_qubit_basis"]

TOL = 1e-9
PAULIS = (X, Y, Z)
I2 = np.eye(2, dtype=complex)

_B = np.array([[1, 0, 0, 1j], [0, 1j, 1, 0], [0, 1j, -1, 0], [1, 0, 0, -1j]]) / math.sqrt(2)
_BD = _B.conj().T
# XX, YY and ZZ are diagonal in the magic basis; their diagonals are ±1.
_DIAG = np.array([np.diag(_BD @ np.kron(p, p) @ _B).real for p in PAULIS])


def canonical(a: float, b: float, c: float) -> np.ndarray:
    """exp(i(a XX + b YY + c ZZ))."""
    return _B @ np.diag(np.exp(1j * (_DIAG.T @ np.array([a, b, c])))) @ _BD


@dataclass
class KAK:
    """``U = e^{i phase} (a1 ⊗ a0) canonical(*point) (b1 ⊗ b0)``.

    ``a1`` and ``b1`` act on the first (most significant) qubit.
    """

    phase: float
    point: tuple[float, float, float]
    a1: np.ndarray
    a0: np.ndarray
    b1: np.ndarray
    b0: np.ndarray

    def matrix(self) -> np.ndarray:
        return (
            np.exp(1j * self.phase)
            * np.kron(self.a1, self.a0)
            @ canonical(*self.point)
            @ np.kron(self.b1, self.b0)
        )


def _factor(k: np.ndarray) -> tuple[np.ndarray, np.ndarray, complex]:
    """Split ``k = phase · (p ⊗ q)`` with ``p`` and ``q`` in SU(2)."""
    m = k.reshape(2, 2, 2, 2).transpose(0, 2, 1, 3).reshape(4, 4)
    u, _, vh = np.linalg.svd(m)
    p = u[:, 0].reshape(2, 2)
    q = vh[0].reshape(2, 2)
    p = p / np.sqrt(np.linalg.det(p))
    q = q / np.sqrt(np.linalg.det(q))
    scale = np.trace(np.kron(p, q).conj().T @ k) / 4
    return p, q, scale


def _orthogonal_eigenvectors(m: np.ndarray) -> np.ndarray:
    """Real orthogonal ``P`` (det 1) diagonalising the symmetric unitary ``m``.

    The real and imaginary parts of ``m`` commute, so a generic real combination
    of them has the common eigenvectors.
    """
    rng = np.random.default_rng(2024)
    for _ in range(100):
        r = rng.uniform()
        _, p = np.linalg.eigh(r * m.real + (1 - r) * m.imag)
        d = p.T @ m @ p
        if np.allclose(d, np.diag(np.diag(d)), atol=1e-10):
            if np.linalg.det(p) < 0:
                p[:, 0] = -p[:, 0]
            return p
    raise np.linalg.LinAlgError("could not diagonalise the KAK matrix")  # pragma: no cover


def kak(u: np.ndarray) -> KAK:
    """The KAK decomposition of a 4×4 unitary, with its point in the Weyl chamber."""
    det = np.linalg.det(u)
    g = np.angle(det) / 4
    up = _BD @ (u * np.exp(-1j * g)) @ _B
    p = _orthogonal_eigenvectors(up.T @ up)
    d = np.sqrt(np.diag(p.T @ up.T @ up @ p))
    k1 = up @ p @ np.diag(d.conj())
    if np.linalg.det(k1).real < 0:
        d[0] = -d[0]
        k1[:, 0] = -k1[:, 0]
    # up = k1 · diag(d) · p.T with k1 and p.T in SO(4).
    theta = np.angle(d)
    offset = theta.mean()
    a, b, c = _DIAG @ (theta - offset) / 4
    left1, left0, s1 = _factor(_B @ k1.real @ _BD)
    right1, right0, s2 = _factor(_B @ p.T @ _BD)
    phase = g + offset + np.angle(s1) + np.angle(s2)
    return _canonicalize(phase, [a, b, c], [left1, left0], [right1, right0])


def _canonicalize(
    phase: float, v: list[float], left: list[np.ndarray], right: list[np.ndarray]
) -> KAK:
    """Move ``v`` into the Weyl chamber, keeping ``U`` unchanged."""
    quarter = math.pi / 4

    def shift(k: int, step: int) -> None:
        # exp(i(v + step·π/2) PP) = exp(i v PP) · (i PP)^step
        nonlocal phase
        v[k] += step * math.pi / 2
        right[0] = PAULIS[k] @ right[0]
        right[1] = PAULIS[k] @ right[1]
        phase -= step * math.pi / 2

    def negate(k1: int, k2: int) -> None:
        # (P ⊗ I) canonical(v) (P ⊗ I) flips the two coordinates P anticommutes with.
        p = PAULIS[3 - k1 - k2]
        v[k1], v[k2] = -v[k1], -v[k2]
        left[0] = left[0] @ p
        right[0] = p @ right[0]

    def swap(k1: int, k2: int) -> None:
        # A quarter turn about the third axis on both qubits exchanges the other two.
        r = _quarter_turn(3 - k1 - k2)
        v[k1], v[k2] = v[k2], v[k1]
        left[0], left[1] = left[0] @ r.conj().T, left[1] @ r.conj().T
        right[0], right[1] = r @ right[0], r @ right[1]

    for k in range(3):
        while v[k] > quarter + TOL:
            shift(k, -1)
        while v[k] <= -quarter + TOL:
            shift(k, 1)
    for _ in range(3):
        for k in range(2):
            if abs(v[k]) < abs(v[k + 1]) - TOL:
                swap(k, k + 1)
    if v[0] < -TOL:
        negate(0, 2)
    if v[1] < -TOL:
        negate(1, 2)
    if abs(v[0] - quarter) < TOL and v[2] < -TOL:
        shift(0, -1)
        negate(0, 2)
    return KAK(phase, (v[0], v[1], v[2]), left[0], left[1], right[0], right[1])


def _quarter_turn(axis: int) -> np.ndarray:
    """exp(-iπ/4 P) for the Pauli on ``axis``."""
    return (I2 - 1j * PAULIS[axis]) / math.sqrt(2)


# -- synthesis --------------------------------------------------------------------

#: A step of a two-qubit sequence, in time order: a native gate
#: (name, params, qubits) or a pending single-qubit unitary ("u", matrix, (qubit,)).
#: Qubit 0 is the first operand of the decomposed gate.
Step = tuple


@dataclass(frozen=True)
class TwoQubitBasis:
    """A native two-qubit gate and how the templates use it."""

    name: str
    #: "cx" for gates equal to CX up to single-qubit gates, "ising" for exp(-iθ/2 P⊗Q).
    kind: str
    #: For "cx": (after1, after0, before1, before0) with
    #: CX ∝ (after1 ⊗ after0) · G · (before1 ⊗ before0).
    locals: tuple[np.ndarray, ...] = ()
    #: For "ising": the Paulis (first, second) of the rotation.
    axes: tuple[int, int] = (2, 2)


_ISING = {"rxx": (0, 0), "ryy": (1, 1), "rzz": (2, 2), "rzx": (2, 0)}


def two_qubit_basis(name: str) -> TwoQubitBasis | None:
    """How to use native gate ``name``, or None when no template uses it yet."""
    if name in _ISING:
        return TwoQubitBasis(name, "ising", axes=_ISING[name])
    try:
        g = gate_matrix(name)
    except (KeyError, ValueError):
        return None
    if g.shape != (4, 4):
        return None
    kg = kak(g)
    if not np.allclose(kg.point, (math.pi / 4, 0, 0), atol=TOL):
        return None
    kc = kak(gate_matrix("cx"))
    # CX ∝ (Ac Ag†) G (Bg† Bc) on each qubit.
    return TwoQubitBasis(
        name,
        "cx",
        locals=(
            kc.a1 @ kg.a1.conj().T,
            kc.a0 @ kg.a0.conj().T,
            kg.b1.conj().T @ kc.b1,
            kg.b0.conj().T @ kc.b0,
        ),
    )


def _cx_templates(point: tuple[float, float, float]) -> list[list[Step]]:
    """Circuits with the fewest CXs whose Weyl point is ``point``.

    Several variants are returned: they differ only by single-qubit rotations
    that commute with a neighbouring CX (Z on a control, X on a target), which
    the outer single-qubit gates absorb, and the cheapest one depends on the
    device's single-qubit basis.
    """
    a, b, c = point
    if np.allclose(point, 0, atol=TOL):
        return [[]]
    if np.allclose(point, (math.pi / 4, 0, 0), atol=TOL):
        return [[("cx", (), (0, 1))]]
    if abs(c) < TOL:
        # CX (Rx(θ) ⊗ Rz(φ)) CX = exp(-iθ/2 XX - iφ/2 ZZ), and either coordinate
        # can take either role. Rx(θ) on the control may be dressed with Z
        # rotations, e.g. into SX Rz(π - θ) SX.
        return [
            [
                ("cx", (), (0, 1)),
                ("u", control, (0,)),
                ("u", rz(-2 * on_target), (1,)),
                ("cx", (), (0, 1)),
            ]
            for on_control, on_target in ((a, b), (b, a))
            for control in (rx(-2 * on_control), _SX @ rz(math.pi - 2 * on_control) @ _SX)
        ]
    half = math.pi / 2
    variants = [
        [
            ("cx", (), (1, 0)),
            ("u", rz(2 * c + half), (0,)),
            ("u", _ry(2 * a + half), (1,)),
            ("cx", (), (0, 1)),
            ("u", _ry(2 * b + half), (1,)),
            ("cx", (), (1, 0)),
        ]
    ]
    # Ry(θ) = Rx(-sπ/2) Rz(sθ) Rx(sπ/2), and the X rotations around the middle
    # CX's target cancel through it.
    for sign in (1, -1):
        variants.append(
            [
                ("cx", (), (1, 0)),
                ("u", rz(2 * c + half), (0,)),
                ("u", rz(sign * (2 * a + half)) @ rx(sign * half), (1,)),
                ("cx", (), (0, 1)),
                ("u", rx(-sign * half) @ rz(sign * (2 * b + half)), (1,)),
                ("cx", (), (1, 0)),
            ]
        )
    return variants


def _ry(theta: float) -> np.ndarray:
    return gate_matrix("ry", (theta,))


def _clifford_map(src: int, dst: int) -> np.ndarray:
    """A single-qubit Clifford V with V P_src V† = P_dst."""
    candidates = [I2, H, S, S.conj().T, S @ H, H @ S, H @ S.conj().T, S.conj().T @ H]
    for v in candidates:
        if np.allclose(v @ PAULIS[src] @ v.conj().T, PAULIS[dst]):
            return v
    raise AssertionError((src, dst))  # pragma: no cover


def _ising_template(point: tuple[float, float, float], basis: TwoQubitBasis) -> list[Step]:
    """exp(i(a XX + b YY + c ZZ)) as one native rotation per non-zero coordinate."""
    steps: list[Step] = []
    first, second = basis.axes
    for axis, coordinate in enumerate(point):
        if abs(coordinate) < TOL:
            continue
        # R_PP(θ) = (V ⊗ W) native(θ) (V ⊗ W)† with V P_first V† = P = W P_second W†.
        v = _clifford_map(first, axis)
        w = _clifford_map(second, axis)
        steps += [
            ("u", v.conj().T, (0,)),
            ("u", w.conj().T, (1,)),
            (basis.name, (-2 * coordinate,), (0, 1)),
            ("u", v, (0,)),
            ("u", w, (1,)),
        ]
    return steps


def _steps_matrix(steps: Sequence[Step]) -> np.ndarray:
    out = np.eye(4, dtype=complex)
    for name, params, qubits in steps:
        if name == "u":
            m = np.kron(params, I2) if qubits == (0,) else np.kron(I2, params)
        else:
            m = gate_matrix(name, params)
            if qubits == (1, 0):
                m = SWAP @ m @ SWAP
        out = m @ out
    return out


def synthesize_2q(
    u: np.ndarray,
    basis: TwoQubitBasis,
    cost_1q: Callable[[np.ndarray], int] | None = None,
) -> list[Step]:
    """``u`` as native ``basis`` gates and single-qubit unitaries, in time order.

    The single-qubit steps are ("u", matrix, (qubit,)) and still need a
    single-qubit synthesis; adjacent ones on the same qubit are merged. When
    ``cost_1q`` gives the gate count of a single-qubit unitary, the template
    variant with the fewest gates in total is chosen.
    """
    target = kak(u)
    if basis.kind == "ising":
        templates = [_ising_template(target.point, basis)]
    else:
        templates = _cx_templates(target.point)
    best: list[Step] | None = None
    best_cost = math.inf
    for template in templates:
        steps = _instantiate(target, template, basis)
        if steps is None:
            continue
        if cost_1q is None:
            return steps
        steps = _dress(steps, cost_1q)
        cost = sum(cost_1q(params) if name == "u" else 1 for name, params, _ in steps)
        if cost < best_cost:
            best, best_cost = steps, cost
    if best is None:
        raise AssertionError(f"no template reaches the point {target.point}")  # pragma: no cover
    return best


def _instantiate(target: KAK, template: list[Step], basis: TwoQubitBasis) -> list[Step] | None:
    """``target`` as ``template`` between single-qubit gates, or None if the points differ."""
    reference = kak(_steps_matrix(template))
    if not np.allclose(reference.point, target.point, atol=1e-7):
        return None
    # u ∝ (A C†) T (D† B) on each qubit, where T is the template.
    steps: list[Step] = [
        ("u", reference.b1.conj().T @ target.b1, (0,)),
        ("u", reference.b0.conj().T @ target.b0, (1,)),
    ]
    for step in template:
        if step[0] == "cx" and basis.name != "cx":
            control, target_qubit = step[2]
            after1, after0, before1, before0 = basis.locals
            steps += [
                ("u", before1, (control,)),
                ("u", before0, (target_qubit,)),
                (basis.name, (), step[2]),
                ("u", after1, (control,)),
                ("u", after0, (target_qubit,)),
            ]
        else:
            steps.append(step)
    steps += [
        ("u", target.a1 @ reference.a1.conj().T, (0,)),
        ("u", target.a0 @ reference.a0.conj().T, (1,)),
    ]
    return _merge_locals(steps)


def _merge_locals(steps: Sequence[Step]) -> list[Step]:
    """Multiply runs of single-qubit unitaries on the same qubit together."""
    out: list[Step] = []
    pending: dict[int, np.ndarray] = {}
    for name, params, qubits in steps:
        if name == "u":
            q = qubits[0]
            pending[q] = params @ pending.get(q, I2)
            continue
        for q in qubits:
            if q in pending:
                out.append(("u", pending.pop(q), (q,)))
        out.append((name, params, qubits))
    for q in sorted(pending):
        out.append(("u", pending[q], (q,)))
    return out


def best_2q_bases(native: Collection[str]) -> list[TwoQubitBasis]:
    """The native two-qubit gates some template can use."""
    return [b for name in sorted(native) if (b := two_qubit_basis(name)) is not None]


def _rotation(axis: int, angle: float) -> np.ndarray:
    return math.cos(angle / 2) * I2 - 1j * math.sin(angle / 2) * PAULIS[axis]


def _oriented(name: str, params: tuple, operands: tuple[int, int]) -> np.ndarray:
    g = gate_matrix(name, params)
    return SWAP @ g @ SWAP if operands == (1, 0) else g


@lru_cache(maxsize=4096)
def _commuting_axes(name: str, params: tuple, operands: tuple[int, int]) -> dict[int, list[int]]:
    """For each operand, the Pauli axes whose rotations commute with the gate."""
    g = _oriented(name, params, operands)
    out: dict[int, list[int]] = {0: [], 1: []}
    for q in (0, 1):
        for axis, p in enumerate(PAULIS):
            m = np.kron(p, I2) if q == 0 else np.kron(I2, p)
            if np.allclose(g @ m, m @ g, atol=1e-12):
                out[q].append(axis)
    return out


def _special_angles(f: Callable[[float], float]) -> list[float]:
    """Angles α where ``f(α) = |m(α)_00|²`` reaches 1, 1/2 or 0.

    For a single-qubit ``m(α)`` that depends on α through one rotation,
    ``f(α) = A + B cos α + C sin α``; those values are where the Euler angle θ
    of ``m(α)`` is 0, π/2 or π, which many bases implement with fewer gates.
    """
    f0, f1, f2 = f(0.0), f(math.pi / 2), f(math.pi)
    a, b, c = (f0 + f2) / 2, (f0 - f2) / 2, f1 - (f0 + f2) / 2
    r = math.hypot(b, c)
    if r < 1e-12:
        return []
    delta = math.atan2(c, b)
    out = []
    for target in (1.0, 0.5, 0.0):
        k = (target - a) / r
        if abs(k) <= 1 + 1e-12:
            x = math.acos(max(-1.0, min(1.0, k)))
            out += [delta + x, delta - x]
    return out


def _candidate_angles(axis: int, before: np.ndarray, after: np.ndarray) -> list[float]:
    """Angles worth trying for a rotation about ``axis`` moved from ``before`` to ``after``."""
    from qevedo.compiler.synthesis.one_qubit import euler_zyz

    out = [0.0]
    # Angles that make θ of either side special.
    out += _special_angles(lambda t: abs((_rotation(axis, t).conj().T @ before)[0, 0]) ** 2)
    out += _special_angles(lambda t: abs((after @ _rotation(axis, t))[0, 0]) ** 2)
    if axis == 2:
        # Z rotations change only the outer Euler angles: cancel one of them.
        # When θ = 0 only the sum of the two Z angles is defined: try it whole.
        _, phi, lam_b = euler_zyz(before)
        _, phi_a, lam = euler_zyz(after)
        bases = (phi, -lam, phi + lam_b, -(phi_a + lam))
        out += [
            base + extra for base in bases for extra in (0.0, math.pi, math.pi / 2, -math.pi / 2)
        ]
    # Rotations are 4π-periodic up to sign, which costs nothing: dedupe modulo 2π.
    unique = {round(math.remainder(angle, 2 * math.pi), 9): angle for angle in out}
    return list(unique.values())


def _dress(steps: list[Step], cost_1q: Callable[[np.ndarray], int]) -> list[Step]:
    """Shift rotations that commute with a two-qubit gate across it to cheapen both sides.

    A rotation about Z on a CX's control (or X on its target) can move from the
    single-qubit gate before the CX to the one after it without changing the
    product; the angle that minimises the two gates' combined cost is kept.
    """
    work = [list(step) for step in steps]
    for _ in range(2):
        i = 0
        while i < len(work):
            name, params, operands = work[i]
            if name == "u":
                i += 1
                continue
            axes = _commuting_axes(name, params, operands)
            for role, q in enumerate(operands):
                for axis in axes[role]:
                    before = _neighbour(work, i, q, -1)
                    if before is None:
                        work.insert(i, ["u", I2, (q,)])
                        before, i = i, i + 1
                    after = _neighbour(work, i, q, 1)
                    if after is None:
                        work.insert(i + 1, ["u", I2, (q,)])
                        after = i + 1
                    b, a = work[before][1], work[after][1]
                    best_cost, best_angle = cost_1q(b) + cost_1q(a), None
                    for angle in _candidate_angles(axis, b, a):
                        r = _rotation(axis, angle)
                        cost = cost_1q(r.conj().T @ b) + cost_1q(a @ r)
                        if cost < best_cost:
                            best_cost, best_angle = cost, angle
                    if best_angle is not None:
                        r = _rotation(axis, best_angle)
                        work[before][1] = r.conj().T @ b
                        work[after][1] = a @ r
            i = _move_paulis(work, i, cost_1q)
            i += 1
    return [tuple(step) for step in work if not (step[0] == "u" and np.allclose(step[1], I2))]


_PAULI_PAIRS = [(p, q) for p in range(4) for q in range(4) if (p, q) != (0, 0)]
_PAULI_1Q = (I2, *PAULIS)
_PAULI_2Q = {(p, q): np.kron(_PAULI_1Q[p], _PAULI_1Q[q]) for p in range(4) for q in range(4)}


@lru_cache(maxsize=4096)
def _pauli_images(name: str, params: tuple, operands: tuple[int, int]) -> tuple:
    """The Pauli products P⊗Q that the gate maps to Pauli products, with their images."""
    g = _oriented(name, params, operands)
    out = []
    for pq in _PAULI_PAIRS:
        image = g @ _PAULI_2Q[pq] @ g.conj().T
        for candidate, matrix in _PAULI_2Q.items():
            if abs(abs(np.vdot(matrix, image)) / 4 - 1) < 1e-9:
                out.append((pq, candidate))
                break
    return tuple(out)


def _move_paulis(work: list, i: int, cost_1q: Callable[[np.ndarray], int]) -> int:
    """Move a Pauli product P⊗Q across Clifford gate ``i`` if that cheapens its neighbours.

    G (P⊗Q) = (P'⊗Q') G when G maps P⊗Q to another Pauli product, as CX maps X⊗I
    to X⊗X. Returns the (possibly shifted) index of the gate.
    """
    name, params, operands = work[i]
    images = _pauli_images(name, tuple(params), tuple(operands))
    if not images:
        return i
    slots = []
    for q in (0, 1):
        before = _neighbour(work, i, q, -1)
        if before is None:
            work.insert(i, ["u", I2, (q,)])
            before, i = i, i + 1
        slots.append(before)
    for q in (0, 1):
        after = _neighbour(work, i, q, 1)
        if after is None:
            work.insert(i + 1, ["u", I2, (q,)])
            after = i + 1
        slots.append(after)
    b0, b1, a0, a1 = (work[j][1] for j in slots)
    best_cost = sum(cost_1q(m) for m in (b0, b1, a0, a1))
    best = None
    for (p, q), (p2, q2) in images:
        moved = (
            _PAULI_1Q[p] @ b0,
            _PAULI_1Q[q] @ b1,
            a0 @ _PAULI_1Q[p2],
            a1 @ _PAULI_1Q[q2],
        )
        cost = sum(cost_1q(m) for m in moved)
        if cost < best_cost:
            best_cost, best = cost, moved
    if best is not None:
        for j, m in zip(slots, best):
            work[j][1] = m
    return i


def _neighbour(steps: list, i: int, q: int, direction: int) -> int | None:
    """The index of the single-qubit step on ``q`` next to step ``i``, if no gate on ``q`` intervenes."""
    j = i + direction
    while 0 <= j < len(steps):
        name, _, operands = steps[j]
        if q in operands:
            return j if name == "u" else None
        j += direction
    return None
