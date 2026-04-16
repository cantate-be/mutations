"""
Tests that are performed:

* `vtjson` is used to check that there is no missing data.
* The case identifiers are consecutive lists that have the correct length.
* The exceptional collections are really exceptional.
* The exceptional collections have the required Gram matrix.
* The exceptional collections are minimal.
* The exceptional collections are very_strong.
* Paths between cases in the relation sections are valid.
* The relations graph for each del Pezzo is a tree.
* The quivers are correct.
* The ranks are correct.
* There is a certificate for each surface.
* Certificates are valid. Either:
  - The symmetry group is trivial or
  - all reflections give the same NCCR or
  - a valid certificate is supplied.
* All Gram matrices are distinct up to rotation.
* All wild card 3 block cases are valid.
"""

import time
from ast import literal_eval
from collections.abc import Sequence
from typing import reveal_type  # noqa: F401
from typing import Any, Literal, TypeAlias, TypedDict, cast

from sage.all import DiGraph, Graph, Matrix  # type: ignore
from vtjson import safe_cast

import mutations as mts

Collection: TypeAlias = list[tuple[int, tuple[int, ...]]]
Certificate: TypeAlias = list[tuple[list[int], list[int]]]
ExtendedCertificate: TypeAlias = (
    Literal["all_reflections_give_same_nccr", "symmetry_group_is_trivial"] | Certificate
)
GramMatrix: TypeAlias = list[list[int]]
Alphas: TypeAlias = list[int]
Ranks: TypeAlias = list[int]
SurfaceName: TypeAlias = Literal[
    "P2", "P1P1", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8"
]

Case: TypeAlias = tuple[Literal[3, 4], int]


class MinimalCollection(TypedDict):
    Mred: GramMatrix
    Qred: list[int]
    alphas: Alphas
    ranks: Ranks
    certificate: None | ExtendedCertificate
    collection: Collection


Relation: TypeAlias = tuple[
    Case,
    Case | Literal["three_block"],
    list[int],
]


class MinimalCollections(TypedDict):
    cases: dict[Case, MinimalCollection]
    relations: list[Relation]


Data: TypeAlias = dict[SurfaceName, MinimalCollections]


def verify_cert(col: mts.ExceptionalCollection, cert: Certificate) -> bool:
    """
    Cert is a sequence of pairs (w, m) where w is a sequence
    of reflections, numbered 0,1,..., ordered as in Sage and
    m is a sequence of left cluster mutations, numbered 1,2,...
    All sequences should be applied from left to right!
    The pairs (w, m) should satisfy w(col) = m(col) where =
    means "having the same NCCR" and the w's together
    should generate the Weyl group.
    """

    def to_elt(seq: Sequence[Any]) -> Any:
        return col.S.weyl_group.from_reduced_word(list(reversed([s + 1 for s in seq])))

    for w, m in cert:
        if not col.apply_reflections(w).has_same_NCCR(col.quiver_mutate(m)):
            return False
        group = col.S.weyl_group.subgroup([to_elt(w) for w, m in cert])
        if len(group) != len(col.S.weyl_group):
            return False
    return True


class TestHarness:
    data: Data

    def __init__(self, minimal_collections_file: str):
        with open(minimal_collections_file) as f:
            txt = f.read()
            self.data = safe_cast(Data, literal_eval(txt))

    def test_all(self) -> bool:
        return all(
            x()
            for x in [
                self.test_keys,
                self.test_anchors,
                self.test_graph_connected,
                self.test_collections,
                self.test_relations,
                self.test_wild_cards,
                self.test_distinct,
            ]
        )

    def test_keys(self) -> bool:
        block3_length = 21
        block4_length = 9
        keys = []
        for surface_name, minimal_collections in self.data.items():
            for k in minimal_collections["cases"]:
                keys.append(k)
        if len(keys) != block3_length + block4_length:
            return False
        keys3 = list(sorted([c for b, c in keys if b == 3]))
        keys4 = list(sorted([c for b, c in keys if b == 4]))
        return keys3 == list(range(1, block3_length + 1)) and keys4 == list(
            range(1, block4_length + 1)
        )

    def test_anchors(self) -> bool:
        for surface, minimal_collections in self.data.items():
            for mc in minimal_collections["cases"].values():
                if mc["certificate"] is not None:
                    break
            else:
                return False
        return True

    def test_graph_connected(self) -> bool:
        for surface, minimal_collections in self.data.items():
            G = DiGraph()
            for case in minimal_collections["cases"]:
                G.add_vertex(case)
                if case[0] == 3:
                    target_3blocks = case
            for relation in minimal_collections["relations"]:
                source = relation[0]
                if isinstance(relation[1], tuple):
                    target = relation[1]
                    G.add_edge(source, target)
                elif relation[1] == "three_block":
                    G.add_edge(source, target_3blocks)
                else:
                    return False
            if not G.is_connected():
                return False
            if not Graph(G).is_tree():
                return False
        return True

    def test_distinct(self) -> bool:
        for surface_name, minimal_collections in self.data.items():
            surface = mts.Surface(surface_name)
            for i, mc in enumerate(minimal_collections["cases"].values()):
                col = mts.ExceptionalCollection(mc["collection"], surface)
                for j, mc1 in enumerate(
                    list(minimal_collections["cases"].values())[i + 1 :]
                ):
                    col1 = mts.ExceptionalCollection(mc1["collection"], surface)
                    if col.same_gram_matrix_up_to_rotation(col1):
                        return False
        return True

    def test_wild_cards(self) -> bool:
        all_collections = {}
        for surface_name, minimal_collections in self.data.items():
            surface = mts.Surface(surface_name)
            for case, mc in minimal_collections["cases"].items():
                col = mts.ExceptionalCollection(mc["collection"], surface)
                all_collections[col] = case

        for surface_name, minimal_collections in self.data.items():
            surface = mts.Surface(surface_name)
            for source_case, target_case, steps in minimal_collections["relations"]:
                if target_case == "three_block":
                    source_collection = minimal_collections["cases"][source_case][
                        "collection"
                    ]
                    source_collection_ = mts.ExceptionalCollection(
                        source_collection, surface
                    )
                    mutated_source_collection = source_collection_.quiver_mutate(steps)
                    block_count = mutated_source_collection.count_blocks()
                    if block_count != 3:
                        return False
                    while True:
                        q = mutated_source_collection.quiver_rank_reducers()
                        if q != []:
                            mutated_source_collection = (
                                mutated_source_collection.quiver_mutate(q[0])
                            )
                        else:
                            break
                    for col in all_collections:
                        if mutated_source_collection.same_gram_matrix_up_to_rotation(
                            col
                        ):
                            break
                    else:
                        return False
        return True

    def test_collections(self) -> bool:
        for surface_name, minimal_collections in self.data.items():
            for mc in minimal_collections["cases"].values():
                if not all(
                    x(surface_name, mc)
                    for x in [
                        self.test_trivial,
                        self.test_collection_exceptional,
                        self.test_collection_very_strong,
                        self.test_collection_gram_matrix,
                        self.test_certificate,
                        self.test_quiver,
                        self.test_ranks,
                        self.test_minimal,
                    ]
                ):
                    return False
        return True

    def test_relations(self) -> bool:
        for surface_name, minimal_collections in self.data.items():
            relations = minimal_collections["relations"]
            if not all(
                self.test_relation(surface_name, relation) for relation in relations
            ):
                return False
        return True

    def test_trivial(self, surface_name: SurfaceName, mc: MinimalCollection) -> bool:
        if not len(mc["alphas"]) == len(mc["Mred"]):
            return False
        return sum(mc["alphas"]) <= 11

    def test_collection_exceptional(
        self, surface_name: SurfaceName, mc: MinimalCollection
    ) -> bool:
        surface = mts.Surface(surface_name)
        col_ = mts.ExceptionalCollection(mc["collection"], surface)
        return col_.is_exceptional()

    def test_collection_very_strong(
        self, surface_name: SurfaceName, mc: MinimalCollection
    ) -> bool:
        surface = mts.Surface(surface_name)
        col_ = mts.ExceptionalCollection(mc["collection"], surface)
        return col_.is_geometric()

    def test_collection_gram_matrix(
        self, surface_name: SurfaceName, mc: MinimalCollection
    ) -> bool:
        surface = mts.Surface(surface_name)
        M = mts.explode(Matrix(mc["Mred"]), mc["alphas"])
        col_ = mts.ExceptionalCollection(mc["collection"], surface)
        gm = col_.gram_matrix()
        return cast(bool, M == gm)

    def test_certificate(
        self, surface_name: SurfaceName, mc: MinimalCollection
    ) -> bool:
        certificate = mc["certificate"]
        surface = mts.Surface(surface_name, sage_numbering=True)
        col = mts.ExceptionalCollection(mc["collection"], surface)
        if certificate is None:
            return True
        elif certificate == "symmetry_group_is_trivial":
            if surface_name not in {"P2", "X1"}:
                return False
        elif certificate == "all_reflections_give_same_nccr":
            for i in range(len(surface.simple_roots)):
                col1 = col.apply_reflection(i)
                if not col.has_same_NCCR(col1):
                    return False
        else:
            return verify_cert(col, certificate)
        return True

    def test_quiver(self, surface_name: SurfaceName, mc: MinimalCollection) -> bool:
        M = mts.explode(Matrix(mc["Mred"]), mc["alphas"])
        Minv = M**-1
        N_ = -Minv + Minv.transpose()
        N = mts.implode(N_, mc["alphas"])
        if len(mc["Qred"]) == 3:
            if [N[0, 1], N[0, 2], N[1, 2]] != mc["Qred"]:
                return False
        elif len(mc["Qred"]) == 6:
            if [N[0, 1], N[0, 2], N[0, 3], N[1, 2], N[1, 3], N[2, 3]] != mc["Qred"]:
                return False
        else:
            return False
        return True

    def test_ranks(self, surface_name: SurfaceName, mc: MinimalCollection) -> bool:
        sums = [sum(mc["alphas"][:i]) for i in range(0, len(mc["alphas"]))]
        for i, s in enumerate(sums):
            if mc["collection"][s][0] != mc["ranks"][i]:
                return False
        return True

    def test_minimal(self, surface_name: SurfaceName, mc: MinimalCollection) -> bool:
        surface = mts.Surface(surface_name)
        col = mts.ExceptionalCollection(mc["collection"], surface)
        return col.quiver_rank_reducers() == []

    def test_relation(self, surface_name: SurfaceName, relation: Relation) -> bool:
        surface = mts.Surface(surface_name)

        source_case = relation[0]
        target_case = relation[1]
        steps = relation[2]

        source_collection = self.data[surface_name]["cases"][source_case]["collection"]
        source_collection_ = mts.ExceptionalCollection(source_collection, surface)
        mutated_source_collection = source_collection_.quiver_mutate(steps)

        if isinstance(target_case, tuple):
            target_collection = self.data[surface_name]["cases"][target_case][
                "collection"
            ]
            target_collection_ = mts.ExceptionalCollection(target_collection, surface)
            if not mutated_source_collection.same_gram_matrix_up_to_rotation(
                target_collection_
            ):
                return False

        elif target_case == "three_block":
            block_count = mutated_source_collection.count_blocks()
            if block_count != 3:
                return False
        else:
            return False
        return True


if __name__ == "__main__":
    t = TestHarness("minimal_collections.txt")
    start = time.time()
    b = t.test_all()
    print(f"Test finished in {time.time() - start:.2f} seconds.")
    if b:
        print("All checks passed")
    else:
        print("Some checks failed")
