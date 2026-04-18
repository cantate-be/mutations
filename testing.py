from ast import literal_eval

import numpy as np
from sage.all import Matrix
from tqdm import tqdm

import mutations as mts
from validate_minimal_collections import TestHarness

# Make sure a run is deterministic...
np.random.seed(0)

t = TestHarness("minimal_collections.txt")

assert t.test_all()

all_matrices = {}

TEST_DATA = "test_data"
TEST_COUNT = 5000

for surface_name, minimal_collections in t.data.items():
    surface = mts.Surface(surface_name)
    for mc in minimal_collections["cases"].values():
        col = mts.ExceptionalCollection(mc["collection"], surface)
        M0 = Matrix(col.gram_matrix(), immutable=True)
        all_matrices[M0] = (surface_name, M0)
        for i in range(1, len(col)):
            M = Matrix(col.rotate_right(i).gram_matrix(), immutable=True)
            all_matrices[M] = (surface_name, M0)


for surface_name in ("P1P1", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8"):
    tqdm.write("Surface name: " + surface_name)
    surface = mts.Surface(surface_name)
    not_found = False
    case_count = len(t.data[surface_name]["cases"])
    matrices = set()
    with open(f"{TEST_DATA}/geometric_collections_len9_{surface_name}.txt") as f:
        lines = np.random.permutation(f.readlines())[0:TEST_COUNT]
        pbar = tqdm(lines)
        for line in pbar:
            l = literal_eval("[" + line + "]")
            col = mts.ExceptionalCollection(l, surface)
            while True:
                while (r := col.quiver_rank_reducer()) is not None:
                    col = col.quiver_mutate(r)
                if (r := col.block_reducer()) is not None:
                    col = col.quiver_mutate(r)
                else:
                    break
            M = Matrix(col.gram_matrix(), immutable=True)
            if M not in all_matrices:
                tqdm.write(str(col) + " was not found")
                not_found = True
                break
            else:
                surface_name_, can_matrix = all_matrices[M]
                assert surface_name_ == surface_name
                matrices.add(can_matrix)
                pbar.set_description(f"Table coverage {len(matrices)}/{case_count}")
    if not_found:
        break
