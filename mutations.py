from __future__ import annotations

import warnings
from collections.abc import Sequence
from fractions import Fraction
from typing import reveal_type  # noqa: F401
from typing import Any, TypeVar, cast

from sage.all import QQ, ZZ, IntegralLattice, Matrix, WeylGroup  # type: ignore

root_system_names = {
    "P2": "A0",
    "P1P1": "A1",
    "X0": "A0",
    "X1": "A0",
    "X2": "A1",
    "X3": "A1xA2",
    "X4": "A4",
    "X5": "D5",
    "X6": "E6",
    "X7": "E7",
    "X8": "E8",
}


def explode(M: Any, alphas: Sequence[int]) -> Any:
    """
    A utility function that explodes M according to the multiplicities
    given by alphas.
    """
    assert M.nrows() == M.ncols() == len(alphas)
    sums = [sum(alphas[:i]) for i in range(1, len(alphas) + 1)]

    def interval(x: int) -> int:
        for i, s in enumerate(sums):
            if s > x:
                return i
        return len(sums) + 1

    size = sum(alphas)
    N = Matrix(QQ, size, size)
    for x in range(size):
        i = interval(x)
        for y in range(size):
            j = interval(y)
            if i == j:
                if x == y:
                    N[x, y] = 1
                else:
                    N[x, y] = 0
            else:
                N[x, y] = M[i, j]
    return N


def implode(M: Any, alphas: Sequence[int]) -> Any:
    """
    The left inverse of explode.
    """
    assert M.nrows() == M.ncols() == sum(alphas)
    sums = [sum(alphas[:i]) for i in range(0, len(alphas))]
    return M[sums, sums]


def gram_matrix(
    ranks: Sequence[int], chis: Sequence[int], require_integral: bool = False
) -> Any:
    """
    A utility function that constructs a Gram matrix out of
    ranks and neighboring chis
    """
    assert len(ranks) == len(chis) + 1
    size = len(ranks)
    M = Matrix(QQ, size, size)
    for i in range(0, size):
        for j in range(0, size):
            if j < i:
                M[i, j] = 0
            elif j == i:
                M[i, j] = 1
            else:
                kappa: Fraction = Fraction(0)
                for k in range(i, j):
                    kappa += Fraction(int(chis[k]), int(ranks[k] * ranks[k + 1]))
                    M[i, j] = int(ranks[i]) * int(ranks[j]) * kappa
                    if require_integral:
                        if not M[i, j].is_integer():
                            raise ValueError("Gram matrix is not integral")
    return M


T = TypeVar("T")


def is_rotation(l1: Sequence[T], l2: Sequence[T]) -> bool:
    """
    Utility function that checks if l2 is obtained from l1 by rotation
    """
    if len(l1) != len(l2):
        return False
    for i in range(0, len(l2)):
        if (l2[i : len(l2)], l2[0:i]) == (
            l1[0 : len(l2) - i],
            l1[len(l2) - i :],
        ):
            return True
    return False


def dot(x: Any, y: Any) -> int:
    return cast(int, x.inner_product(y))


class Surface:
    __instances: dict[tuple[str, bool], Surface] = {}

    __name__: str
    K0rk: int
    Pic: Any
    K: Any
    Kinv: Any
    Ksquare: int
    simple_roots: list[Any]
    weyl_group: Any
    O: K0Element
    o: K0Element

    def __new__(cls, name: str, sage_numbering: bool = False) -> Surface:
        if (name, sage_numbering) in cls.__instances:
            return cls.__instances[(name, sage_numbering)]
        else:
            return super().__new__(cls)

    def __init__(self, name: str, sage_numbering: bool = False) -> None:
        if (name, sage_numbering) in self.__instances:
            return
        self.__name__ = name
        self.sage_numbering = sage_numbering
        self.weyl_group = None
        if name == "P1P1":
            self.__init_P1P1()
        elif name in ["X0", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8"]:
            self.__init_X(int(name[1]))
        elif name == "P2":
            self.__init_X(0)
        else:
            raise ValueError(f"Invalid surface {name}")

        self.O = K0Element(1, self.Pic(0), 1, self)
        self.o = K0Element(0, self.Pic(0), 1, self)
        self.__instances[(name, sage_numbering)] = self

    @classmethod
    def _clear(cls) -> None:
        cls.__instances = {}

    def __init_P1P1(self) -> None:
        self.K0rk = 4
        gm = Matrix(ZZ, self.K0rk - 2, self.K0rk - 2)
        gm[1, 0] = 1
        gm[0, 1] = 1
        gram = Matrix(ZZ, gm)
        self.Pic = IntegralLattice(gram)

        self.K = self.Pic((-2, -2))
        self.Kinv = -self.K
        self.Ksquare = dot(self.K, self.K)
        assert self.Ksquare == 12 - self.K0rk
        self.simple_roots = [self.Pic((1, -1))]
        if self.sage_numbering:
            self.weyl_group = WeylGroup(
                root_system_names[self.__name__], implementation="permutation"
            )

    def __init_X(self, n: int) -> None:
        self.K0rk = n + 3
        gm = Matrix(ZZ, self.K0rk - 2, self.K0rk - 2)
        gm[0, 0] = 1
        for i in range(1, self.K0rk - 2):
            gm[i, i] = -1
        gram = Matrix(ZZ, gm)
        self.Pic = IntegralLattice(gram)

        self.K = self.Pic((-3,) + n * (1,))
        self.Kinv = -self.K
        self.Ksquare = dot(self.K, self.K)
        H = self.Pic((1,) + n * (0,))
        E: dict[int, Any] = {}
        for i in range(1, n + 1):
            E[i] = (n + 1) * [0]
            E[i][i] = 1
            E[i] = self.Pic(E[i])
        if n in {0, 1}:
            self.simple_roots = []
        elif n == 2:
            self.simple_roots = [E[1] - E[2]]
        elif n >= 3:
            self.simple_roots = [H - E[1] - E[2] - E[3]]
            for i in range(2, n + 1):
                self.simple_roots.append(E[i - 1] - E[i])
            if self.sage_numbering:
                sr = self.simple_roots
                if n == 4:
                    sr[1], sr[3] = sr[3], sr[1]
                elif n == 5:
                    sr[0], sr[1], sr[2], sr[3] = sr[1], sr[2], sr[3], sr[0]
                elif n in (6, 7, 8):
                    sr[0], sr[1] = sr[1], sr[0]
                self.weyl_group = WeylGroup(
                    root_system_names["X" + str(n)], implementation="permutation"
                )
            assert len(self.simple_roots) == n

    def __str__(self) -> str:
        return self.__name__


class K0Element:
    rank: int
    c1: Any
    _chi: int
    surface: Surface

    def __init__(self, rank: int, c1: Any, chi: int, surface: Surface):
        self.rank = rank
        self._chi = chi
        self.S = surface
        self.c1 = self.S.Pic(c1) if isinstance(c1, tuple) else c1
        self.degree = dot(self.S.Kinv, self.c1)

    def chi(self, other: K0Element) -> int:
        """
        Computes chi(self, other)
        """
        return (
            self.rank * other._chi
            + other.rank * self._chi
            - self.rank * other.rank
            - dot(self.c1, other.c1)
            - other.rank * self.degree
        )

    def __str__(self) -> str:
        return str((self.rank, tuple(self.c1), self._chi))


class Quiver:

    N: Any

    def __init__(self, M: Any):
        Minv = M**-1
        self.N = -Minv + Minv.transpose()

    def __str__(self) -> str:
        return str(self.N)

    def has_zero(self) -> bool:
        N = self.N
        size = N.nrows()
        for i in range(0, size):
            for j in range(i + 1, size):
                if abs(N[i, j]) == 0:
                    return True
        return False

    def blocks(self) -> list[list[int]]:
        N = self.N
        size = N.nrows()
        blocks_: list[list[int]] = []
        for i in range(0, size):
            found = False
            for b in blocks_:
                for j in b:
                    if N[i, j] == 0:
                        b.append(i)
                        found = True
                        break
                if found:
                    break
            if not found:
                blocks_.append([i])
        return blocks_

    def is_complete(self) -> bool:
        return not self.has_zero()

    def is_block_complete(self) -> bool:
        N = self.N
        blocks_ = self.blocks()
        for i, b in enumerate(blocks_):
            for c in blocks_[i:]:
                mult = N[b[0], c[0]]
                for bb in b:
                    for cc in c:
                        if b != c:
                            assert N[bb, cc] != 0
                        if b is c:
                            assert N[bb, cc] == 0
                        if N[bb, cc] != mult:
                            return False
        return True


class ExceptionalObject:
    rank: int
    c1: Any
    __chi: None | int
    __degree: None | int
    __slope: None | Fraction
    __hash: None | int
    surface: Surface

    def __init__(self, rank: int, c1: Any, surface: Surface, chi: None | int = None):
        self.rank = rank
        self.S = surface
        if isinstance(c1, tuple) and len(c1) != self.S.K0rk - 2:
            raise ValueError("c1 argument has wrong length")
        self.c1 = self.S.Pic(c1) if isinstance(c1, tuple) else c1
        self.__degree = None
        self.__slope = None
        self.__hash = None
        if self.rank == 0:
            if chi is None:
                raise ValueError(
                    "For rank zero exceptional objects you have to specify chi"
                )
            else:
                self.__chi = chi
        else:
            if chi is not None:
                warnings.warn("Supplied chi value ignored.")
            self.__chi = None

    @property
    def _chi(self) -> int:
        if self.__chi is not None:
            return self.__chi
        else:
            defect = (
                -self.rank * self.rank
                - dot(self.c1, self.c1)
                - self.rank * self.degree
                - 1
            )
            chi, r = divmod(defect, 2 * self.rank)
            if r != 0:
                raise ValueError("self is not an exceptional object")
            return -chi

    @property
    def degree(self) -> int:
        if self.__degree is None:
            self.__degree = dot(self.S.Kinv, self.c1)
        return self.__degree

    @property
    def slope(self) -> Fraction:
        if self.__slope is None:
            self.__slope = Fraction(int(self.degree), int(self.rank))
        return self.__slope

    def is_exceptional(self) -> bool:
        """
        Checks if an object is _numerically_ exceptional.
        """
        if self.rank == 0:
            return -dot(self.c1, self.c1) == 1
        elif self.__chi is None:
            try:
                self._chi  # tries to compute self.__chi and saves it
                return True
            except ValueError:
                return False
        else:
            return (
                2 * self.rank * self._chi
                - self.rank * self.rank
                - dot(self.c1, self.c1)
                + self.rank * self.degree
                == 1
            )

    def is_pair(self, other: ExceptionalObject) -> bool:
        return other.chi(self, assume_pair=False) == 0

    def k0_element(self) -> K0Element:
        return K0Element(self.rank, self.c1, self._chi, self.S)

    def chi(self, other: ExceptionalObject, assume_pair: bool = False) -> int:
        """
        Computes chi(self, other) under the assumption that
        (self, other) is an exceptional pair when assume_pair is True.
        """
        if assume_pair:
            diff = self.rank * other.c1 - other.rank * self.c1
            return dot(self.S.Kinv, diff)
        else:
            k0_self = self.k0_element()
            k0_other = other.k0_element()
            return k0_self.chi(k0_other)

    def derived_shift(self, n: int) -> ExceptionalObject:
        if n % 2 == 0:
            return self
        else:
            if self.rank != 0:
                return self.__class__(-self.rank, -self.c1, self.S)
            else:
                return self.__class__(0, -self.c1, self.S, -self._chi)

    def left_mutate(self, E: ExceptionalObject) -> ExceptionalObject:
        """
        Computes L_E(self)
        """
        ch = E.chi(self, assume_pair=True)  # chi(E,self)
        new_chi: None | int = None
        new_rank = self.rank - ch * E.rank
        if new_rank == 0:
            new_chi = self._chi - ch * E._chi
        return self.__class__(new_rank, self.c1 - ch * E.c1, self.S, new_chi)

    def right_mutate(self, F: ExceptionalObject) -> ExceptionalObject:
        """
        Computes R_F(self)
        """
        ch = self.chi(F, assume_pair=True)  # chi(self,F)
        new_chi: None | int = None
        new_rank = self.rank - ch * F.rank
        if new_rank == 0:
            new_chi = self._chi - ch * F._chi
        return self.__class__(new_rank, self.c1 - ch * F.c1, self.S, new_chi)

    def tensor_by_line_bundle(self, L: Any) -> ExceptionalObject:
        if self.rank != 0:
            return self.__class__(self.rank, self.c1 + self.rank * L, self.S)
        else:
            return self.__class__(0, self.c1, self.S, self._chi + dot(self.c1, L))

    def twist_by_K(self, n: int = 1) -> ExceptionalObject:
        return self.tensor_by_line_bundle(n * self.S.K)

    def twist_by_Kinv(self, n: int = 1) -> ExceptionalObject:
        return self.twist_by_K(-n)

    def line_bundle_difference(self, other: ExceptionalObject) -> Any:
        """
        Finds the line bundle L such that L otimes self = other
        """
        if self.rank != other.rank:
            return None

        if self.rank == 0:
            raise ValueError("This function does not work for rank zero objects")

        nL = (other.c1 - self.c1) / self.rank
        if nL not in self.S.Pic:
            return None
        return self.S.Pic(nL)  # the cast speeds up things enormously

    def make_rank_positive(self) -> ExceptionalObject:
        if self.rank >= 0:
            return self
        else:
            return self.__class__(-self.rank, -self.c1, self.S)

    def apply_reflection(self, i: int) -> ExceptionalObject:
        s = self.S.simple_roots[i]
        if self.rank != 0:
            return self.__class__(self.rank, self.c1 + dot(self.c1, s) * s, self.S)
        else:
            return self.__class__(0, self.c1 + dot(self.c1, s) * s, self.S, self._chi)

    def apply_reflections(self, r: list[int]) -> ExceptionalObject:
        obj = self
        for j in r:
            obj = obj.apply_reflection(j)
        return obj

    def __eq__(self, other: object, /) -> bool:
        if not isinstance(other, self.__class__):
            return False
        if self.rank != 0:
            return (self.S, self.rank, self.c1) == (other.S, other.rank, other.c1)
        else:
            return (self.S, self.c1, self._chi) == (other.S, other.c1, other._chi)

    def __hash__(self) -> int:
        if self.__hash is not None:
            return self.__hash
        if self.rank != 0:
            self.__hash = hash((self.S, self.rank, tuple(self.c1)))
        else:
            self.__hash = hash((self.S, tuple(self.c1), self._chi))
        return self.__hash

    def __str__(self) -> str:
        try:
            return str((self.rank, tuple(self.c1), self._chi))
        except ValueError:
            return str((self.rank, tuple(self.c1), "invalid_chi"))

    __repr__ = __str__


class ExceptionalCollection:

    objects: Sequence[ExceptionalObject]
    S: Surface
    __hash: None | int

    def __init__(
        self,
        objects: Sequence[
            ExceptionalObject
            | tuple[int, tuple[int, ...]]
            | tuple[int, tuple[int, ...], int]
        ],
        surface: Surface,
    ):
        self.S = surface
        self.objects = []
        self.__hash = None
        for o in objects:
            if isinstance(o, ExceptionalObject):
                if o.S != self.S:
                    raise ValueError(f"{o} has an incorrect surface")
                self.objects.append(o)
            elif len(o) == 2:
                E = ExceptionalObject(o[0], o[1], self.S)
                self.objects.append(E)
            else:
                E = ExceptionalObject(o[0], o[1], self.S, o[2])
                self.objects.append(E)
        self.objects = tuple(self.objects)

    def is_exceptional(self) -> bool:
        exceptional = True
        N = len(self)
        for i in range(0, N):
            if not self.objects[i].is_exceptional():
                return False
            for j in range(i + 1, N):
                if not self.objects[i].is_pair(self.objects[j]):
                    exceptional = False
                    break
            if not exceptional:
                break
        return exceptional

    def tensor_by_line_bundle(self, L: Any) -> ExceptionalCollection:
        return self.__class__(
            tuple(o.tensor_by_line_bundle(L) for o in self.objects), self.S
        )

    def make_ranks_positive(self) -> ExceptionalCollection:
        if all([o.rank >= 0 for o in self.objects]):
            return self

        return self.__class__(
            tuple(o.make_rank_positive() for o in self.objects), self.S
        )

    def count_blocks(self) -> int:
        s = self.slopes()
        s += [s[0] + int(self.S.Ksquare)]
        return sum([s[i] != s[i + 1] for i in range(0, len(self))])

    def blocks(self) -> list[list[int]]:
        N = len(self)
        helix = Helix(self)
        current = 0
        if helix[current - 1].slope == helix[current].slope:
            current -= 1

        blocks = []
        current_block = []
        for i in range(current, current + len(self)):
            current_block.append(i % N)
            if helix[i].slope == helix[i + 1].slope:
                continue
            blocks.append(current_block)
            current_block = []
        return blocks

    def block_reducer(self) -> None | int:
        """
        Returns i such that quiver_mutate(i)
        reduces the length of the shortest long edge
        which occurs in a parallel pair. Returns None
        if self is block complete.
        """
        me = self.make_ranks_positive()
        base_block_count = me.count_blocks()
        base_rank = me.rank_sum()
        blocks = me.blocks()
        for b in blocks:
            if len(b) == 1:
                m = b[0] if b[0] > 0 else len(me)
                new_col = me.quiver_mutate(m)
                if (
                    new_col.count_blocks() < base_block_count
                    and new_col.rank_sum() <= base_rank
                ):
                    return m
        min_block_length = None
        min_block = None
        for b in blocks:
            if len(b) == 1:
                continue
            if min_block_length is not None and len(b) >= min_block_length:
                continue
            m = b[0] if b[0] > 0 else len(me)
            new_col = me.quiver_mutate(m)
            if new_col.rank_sum() > base_rank:
                continue
            if new_col.count_blocks() > base_block_count:
                continue
            min_block_length = len(b)
            min_block = b
        if min_block is not None:
            return min_block[0] if min_block[0] > 0 else len(me)
        return None

    def gram_matrix(self) -> Any:
        N = len(self)
        M = Matrix(ZZ, N, N)
        for i in range(0, N):
            for j in range(i, N):
                if i == j:
                    M[i, j] = 1
                else:
                    M[i, j] = self.objects[i].chi(self.objects[j], assume_pair=True)
        return M

    def quiver(self) -> Quiver:
        return Quiver(self.gram_matrix())

    def apply_reflection(self, i: int) -> ExceptionalCollection:
        return self.__class__(
            tuple(o.apply_reflection(i) for o in self.objects), self.S
        )

    def apply_reflections(self, r: list[int]) -> ExceptionalCollection:
        return self.__class__(
            tuple(o.apply_reflections(r) for o in self.objects), self.S
        )

    def slopes(self) -> list[Fraction]:
        return [o.slope for o in self.objects]

    def is_geometric(self) -> bool:
        N = len(self)
        if any(obj.rank <= 0 for obj in self.objects):
            return False
        mu = self.slopes()
        if any(mu[i] > mu[i + 1] for i in range(N - 1)):
            return False
        return mu[-1] - mu[0] <= int(self.S.Ksquare)

    def rotate_left(self, n: None | int = None) -> ExceptionalCollection:
        E5_tw = self.objects[-1].twist_by_K()  # uses that dim S is even
        if n is None:
            return self.__class__((E5_tw,) + tuple(self.objects[:-1]), self.S)
        else:
            col = self
            for _ in range(n):
                col = col.rotate_left()
            return col

    def rotate_right(self, n: None | int = None) -> ExceptionalCollection:
        E1_tw = self.objects[0].twist_by_Kinv()  # uses that dim S is even
        if n is None:
            return self.__class__(tuple(self.objects[1:]) + (E1_tw,), self.S)
        else:
            col = self
            for _ in range(n):
                col = col.rotate_right()
            return col

    def mutate(self, i: int) -> ExceptionalCollection:
        N = len(self)
        if not -N <= i <= N or i == 0:
            raise ValueError(f"length={N}, so mutating at {i} is not valid")
        col = self
        idx = abs(i) - 1
        if idx == N - 1:
            col = col.rotate_right()
            idx = N - 2
        E, F = col.objects[idx], col.objects[idx + 1]
        if i > 0:
            Fprime = F.left_mutate(E)
            new_objects = list(col.objects)
            new_objects[idx] = Fprime
            new_objects[idx + 1] = E  # make tuple
        else:
            Eprime = E.right_mutate(F)
            new_objects = list(col.objects)
            new_objects[idx] = F
            new_objects[idx + 1] = Eprime  # make tuple

        col = self.__class__(tuple(new_objects), self.S)
        if abs(i) == N:
            col = col.rotate_left()

        return col

    def is_fixed_by(self, i: int) -> bool:
        """
        Checks if sigma_i does not change self
        """
        i = abs(i) - 1
        return self[i].chi(Helix(self)[i + 1]) == 0

    def quiver_mutate_left(self, i: int) -> ExceptionalCollection:
        """
        Uses indexing from 1 to len(self)
        """
        N = len(self)
        if not 1 <= i <= N:
            raise ValueError(f"length={N} so left cluster mutating at {i} is invalid")
        col = self.mutate(i)
        while col.is_fixed_by(i):
            i -= 1
            if i == 0:
                i = N
            col = col.mutate(i)
        while True:
            col = col.make_ranks_positive()
            if col.is_geometric():
                return col
            i -= 1
            if i == 0:
                i = N
            col = col.mutate(i)

    def quiver_mutate_right(self, i: int) -> ExceptionalCollection:
        """
        Uses indexing from 1 to len(self)
        """
        N = len(self)
        if not 1 <= i <= N:
            raise ValueError(f"length={N} so right cluster mutating at {i} is invalid")

        col = self.mutate(-i)
        while col.is_fixed_by(i):
            i += 1
            if i == N + 1:
                i = 1
            col = col.mutate(-i)
        while True:
            col = col.make_ranks_positive()
            if col.is_geometric():
                return col
            i += 1
            if i == N + 1:
                i = 1
            col = col.mutate(-i)

    def quiver_mutate(self, i: int | Sequence[int]) -> ExceptionalCollection:
        if isinstance(i, Sequence):
            col = self
            for ii in i:
                col = col.quiver_mutate(ii)
            return col
        if i > 0:
            return self.quiver_mutate_left(i)
        else:
            return self.quiver_mutate_right(-i)

    def quiver_rank_reducer(self) -> None | int:
        N = len(self)
        me = self.make_ranks_positive()
        base_sum = me.rank_sum()
        for i in range(1, N + 1):
            new = me.quiver_mutate(i)
            if new.rank_sum() < base_sum:
                return i
        return None

    def quiver_rank_reducers(self) -> list[int]:
        N = len(self)
        me = self.make_ranks_positive()
        base_sum = me.rank_sum()
        vals: list[int] = []

        for i in range(1, N + 1):
            new = me.quiver_mutate(i)
            if new.rank_sum() < base_sum:
                vals.append(i)
        return vals

    def same_gram_matrix_up_to_rotation(self, target: ExceptionalCollection) -> bool:
        gm_target = target.gram_matrix()
        for i in range(0, len(target)):
            gm_mutated = self.rotate_right(i).gram_matrix()
            if gm_target == gm_mutated:
                return True
        else:
            return False

    def rank_sum(self) -> int:
        return sum(self.ranks())

    def ranks(self) -> list[int]:
        return [o.rank for o in self.objects]

    def normalize_slopes(self) -> ExceptionalCollection:
        """
        Rotates self in such a way that self[0] is the first object
        in the helix with slope >=0. Also make ranks positive.
        """
        col = self.make_ranks_positive()
        if col.objects[0].slope < 0:
            while col.objects[0].slope < 0:
                col = col.rotate_right()
        else:
            while col.objects[0].slope >= 0:
                previous_col = col
                col = col.rotate_left()
            col = previous_col
        return col

    def has_same_NCCR(self, other: ExceptionalCollection) -> bool:
        me = self.make_ranks_positive()
        other = other.make_ranks_positive()

        if me.rank_sum() != other.rank_sum():
            return False
        if not is_rotation(me.ranks(), other.ranks()):
            return False
        ref = self[0]
        other_: None | set[ExceptionalObject] = None
        for o in other.objects:
            L = ref.line_bundle_difference(o)
            if L is None:
                continue
            selfL = me.tensor_by_line_bundle(L)
            # Checks if me and selfL agree up to rotation
            # and permuting blocks
            if other_ is None:
                other_ = set(other.normalize_slopes().objects)
            if set(selfL.normalize_slopes().objects) == other_:
                return True
        return False

    def __getitem__(self, n: int, /) -> ExceptionalObject:
        return self.objects[n]

    def __hash__(self) -> int:
        if self.__hash is not None:
            return self.__hash
        self.__hash = hash((self.S, self.objects))
        return self.__hash

    def __eq__(self, other: object, /) -> bool:
        if not isinstance(other, self.__class__):
            return False
        return (self.S, self.objects) == (self.S, other.objects)

    def __len__(self) -> int:
        return len(self.objects)

    def __str__(self) -> str:
        return "[ " + ", ".join(map(str, self.objects)) + " ]"

    __repr__ = __str__


class Helix:

    col: ExceptionalCollection
    __computed_cols: dict[int, ExceptionalObject]

    def __init__(self, col: ExceptionalCollection):
        self.col = col
        self.__computed_cols = {}

    def __getitem__(self, i: int, /) -> ExceptionalObject:
        """
        Returns the object with index n in the helix generated
        by self.
        """
        if i in self.__computed_cols:
            return self.__computed_cols[i]
        N = len(self.col)
        q, r = divmod(i, N)
        ret = self.col.objects[r].twist_by_Kinv(q)
        self.__computed_cols[i] = ret
        return ret
