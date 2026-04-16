from ast import literal_eval

from sage.all import Matrix
from tqdm import tqdm

import mutations as mts
from validate_minimal_collections import TestHarness

t = TestHarness("minimal_collections.txt")

assert t.test_all()

all_matrices = set()

TEST_DATA = "test_data"

for surface_name, minimal_collections in t.data.items():
    surface = mts.Surface(surface_name)
    for mc in minimal_collections["cases"].values():
        col = mts.ExceptionalCollection(mc["collection"], surface)
        for i in range(len(col)):
            M = Matrix(col.rotate_right(i).gram_matrix(), immutable=True)
            all_matrices.add(M)


for surface_name in ("P1P1", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8"):
    tqdm.write("Surface name: " + surface_name)
    surface = mts.Surface(surface_name)
    not_found = False
    count = 0
    with open(f"{TEST_DATA}/geometric_collections_len9_{surface_name}.txt") as f:
        lines = f.readlines()
        for line in (pbar := tqdm(lines)):

            l = literal_eval("[" + line + "]")
            col = mts.ExceptionalCollection(l, surface)
            if not col.quiver().is_block_complete():
                continue
            if col.quiver_rank_reducers() != []:
                continue
            count += 1
            pbar.set_description("Minimal: " + str(count))
            if Matrix(col.gram_matrix(), immutable=True) not in all_matrices:
                tqdm.write(str(col) + " was not found")
                not_found = True
                break
    if not_found:
        break
