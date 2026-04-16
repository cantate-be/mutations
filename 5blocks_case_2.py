import ast
import math
import os
import time
from itertools import product

import numpy as np
from sage.all import Integer
from tqdm import tqdm

from mutations import explode, gram_matrix
from sum_lemma import sum_lemma_max

# Make sure a run is deterministic...
np.random.seed(0)

DATA_DIR = "5blocks_case_2_data"

if not os.path.exists(DATA_DIR):
    os.mkdir(DATA_DIR)


def is_unimodal4(l):
    """
    A helper that checks if r1,r2,r3,r4 satifies condition (13.65).
    Very inefficient implementation... Note also that "unimodal"
    is a misnomer.
    """
    r1, r2, r3, r4 = l
    if (
        r1 == r2 == r3 == r4
        or r1 > r2 == r3 == r4
        or r1 > r2 > r3 == r4
        or r1 > r2 > r3 > r4
        or r1 == r2 == r3 < r4
        or r1 > r2 == r3 < r4
        or r1 > r2 > r3 < r4
        or r1 == r2 < r3 < r4
        or r1 > r2 < r3 < r4
        or r1 < r2 < r3 < r4
    ):
        return True
    else:
        return False


def is_unimodal3(l):
    r1, r2, r3 = l
    return r1 >= r2 <= r3


# compute the list of possible alpha0,...,alpha4, K2, chi23
# and write it to a file
def valid_settings():
    """ """
    with open(f"{DATA_DIR}/5blocks_case_2_stage_0.txt", "w") as f:
        for chi23 in range(1, 9):
            for K2 in range(1, 10):
                for alpha1, alpha2 in [(1, 1), (1, 2), (1, 3), (2, 1), (3, 1)]:
                    for alpha3, alpha4 in [(1, 1), (1, 2), (1, 3), (2, 1), (3, 1)]:
                        alpha0 = 12 - K2 - alpha1 - alpha2 - alpha3 - alpha4
                        if alpha0 <= 0:
                            continue
                        out = {
                            "chi23": chi23,
                            "K2": K2,
                            "alpha1": alpha1,
                            "alpha2": alpha2,
                            "alpha3": alpha3,
                            "alpha4": alpha4,
                            "alpha0": alpha0,
                        }
                        f.write(str(out) + "\n")


# Finds an upper bound for r0 using the iteration procedure
# described in the ms.
def process_setting_stage_1(s):
    K2 = Integer(s["K2"])
    chi23 = Integer(s["chi23"])
    alpha0 = Integer(s["alpha0"])
    alpha1 = Integer(s["alpha1"])
    alpha2 = Integer(s["alpha2"])
    alpha3 = Integer(s["alpha3"])
    alpha4 = Integer(s["alpha4"])
    assert 12 == K2 + alpha0 + alpha1 + alpha2 + alpha3 + alpha4

    def f1(x):
        return 1 / Integer(x)

    def f2(x):
        return chi23 / Integer(x)

    g = (f1, f2, f1)
    m = sum_lemma_max(g, K2)
    gap = K2 - m
    r0_upperbound_best = Integer(
        math.floor(math.sqrt((16 / (gap * alpha0)) + 1))
    )  # +1 because of possible rounding errors
    while True:
        m = -1
        for r1, r2, r3 in tqdm(
            list(product(range(1, r0_upperbound_best + 1), repeat=3)),
            smoothing=0.0,
            leave=False,
        ):
            if not is_unimodal3([r1, r2, r3]):
                continue
            v = Integer(1) / (r1 * r2) + Integer(chi23) / (r2 * r3)
            if v >= K2:
                continue
            gap = K2 - v
            r4_lower_bound = (1 / (gap * r3)).floor() + 1
            u = (
                Integer(1) / (r1 * r2)
                + Integer(chi23) / (r2 * r3)
                + Integer(1) / (r3 * r4_lower_bound)
            )
            assert u <= K2
            if u > m:
                m = u
        gap = K2 - m
        r0_upperbound_new_try = Integer(
            math.floor(math.sqrt((16 / (gap * alpha0)) + 1))
        )  # +1 because of possible rounding errors
        if r0_upperbound_new_try >= r0_upperbound_best:
            break
        else:
            r0_upperbound_best = r0_upperbound_new_try
    return r0_upperbound_best


def process_settings_stage_1():
    with open(f"{DATA_DIR}/5blocks_case_2_stage_0.txt") as f:
        lines = f.readlines()
        # To make uniform progress...
        lines = np.random.permutation(lines)
    with open(f"{DATA_DIR}/5blocks_case_2_stage_1.txt", "w") as f:
        for line in tqdm(lines, smoothing=0.0):
            s = ast.literal_eval(line)
            r0_upperbound = process_setting_stage_1(s)
            s["r0_upperbound"] = r0_upperbound
            f.write(str(s) + "\n")


# Now find those r1,r2,r3,r4 for which the quiver
# has no backward arrows
def process_setting_stage_2(s):
    sols = []
    K2 = Integer(s["K2"])
    chi23 = Integer(s["chi23"])
    alpha0 = Integer(s["alpha0"])
    alpha1 = Integer(s["alpha1"])
    alpha2 = Integer(s["alpha2"])
    alpha3 = Integer(s["alpha3"])
    alpha4 = Integer(s["alpha4"])
    r0_upperbound = Integer(s["r0_upperbound"])
    assert 12 == K2 + alpha0 + alpha1 + alpha2 + alpha3 + alpha4
    p = product(range(1, r0_upperbound + 1), repeat=2)
    for r1, r2 in tqdm(list(p), leave=False, smoothing=0.0):
        s2 = Integer(1) / (r1 * r2)
        for r3 in range(1, r0_upperbound + 1):
            s3 = s2 + Integer(chi23) / (r2 * r3)
            if s3 >= K2:
                continue
            ranks_ = (r1, r2, r3)
            chis_ = (1, chi23)
            try:
                gram_matrix(ranks_, chis_, require_integral=True)
            except ValueError:
                continue
            for r4 in range(1, r0_upperbound + 1):
                s4 = s3 + Integer(1) / (r3 * r4)
                if s4 < K2:
                    ranks = (r1, r2, r3, r4)
                    if not is_unimodal4(ranks):
                        continue
                    chis = (1, chi23, 1)
                    try:
                        M = gram_matrix(ranks, chis, require_integral=True)
                    except ValueError:
                        continue
                    assert M.nrows() == 4
                    alphas = [alpha1, alpha2, alpha3, alpha4]
                    alpha_sum = sum(alphas)
                    M = explode(M, alphas)
                    N = M**-1
                    good = True
                    for i in range(alpha_sum):
                        for j in range(i + 1, alpha_sum):
                            if N[i, j] > 0:  # positive value means relation from i->j
                                good = False
                                break
                        if not good:
                            break
                    if good:
                        sols.append([r1, r2, r3, r4])
    return sols


def process_settings_stage_2():
    with open(f"{DATA_DIR}/5blocks_case_2_stage_1.txt") as f:
        lines = f.readlines()
        # To make uniform progress...
        lines = np.random.permutation(lines)
    with open(f"{DATA_DIR}/5blocks_case_2_stage_2.txt", "w") as f:
        for line in tqdm(lines, smoothing=0.0):
            s = ast.literal_eval(line)
            sols = process_setting_stage_2(s)
            s["sols"] = sols
            f.write(str(s) + "\n")
            f.flush()


# now we finally study the full Gram matrix
def process_setting_stage_3(s):
    sols = []
    K2 = Integer(s["K2"])
    chi23 = Integer(s["chi23"])
    alpha0 = Integer(s["alpha0"])
    alpha1 = Integer(s["alpha1"])
    alpha2 = Integer(s["alpha2"])
    alpha3 = Integer(s["alpha3"])
    alpha4 = Integer(s["alpha4"])
    assert 12 == K2 + alpha0 + alpha1 + alpha2 + alpha3 + alpha4

    r0_upperbound = Integer(s["r0_upperbound"])

    def to_Integer(x):
        if isinstance(x, (tuple, list)):
            return [to_Integer(y) for y in x]
        else:
            return Integer(x)

    sols = to_Integer(s["sols"])

    data = []

    for r1, r2, r3, r4 in tqdm(sols, smoothing=0.0, leave=False):
        assert isinstance(r1, Integer)
        for r0 in range(max(r1, r2, r3, r3), r0_upperbound + 1):
            gap = (
                K2
                - Integer(1) / (r1 * r2)
                - Integer(chi23) / (r2 * r3)
                - Integer(1) / (r3 * r4)
            )
            chi40_upperbound1 = (gap * r4 * r0).floor()
            chi40_upperbound2 = (
                (Integer(8) * r4) / (alpha0 * r0)
            ).floor()  # (13.66) in the manuscript
            chi40_upperbound = min(chi40_upperbound1, chi40_upperbound2)
            assert chi40_upperbound / (r4 * r0) <= gap
            for chi40 in range(1, chi40_upperbound + 1):
                chi40 = Integer(chi40)
                gap2 = gap - chi40 / (r4 * r0)
                chi01 = gap2 * (r0 * r1)
                assert (
                    Integer(1) / (r1 * r2)
                    + Integer(chi23) / (r2 * r3)
                    + Integer(1) / (r3 * r4)
                    + chi40 / (r4 * r0)
                    + chi01 / (r0 * r1)
                    == K2
                )
                if not chi01.is_integer():
                    continue

                try:
                    M = gram_matrix(
                        (r0, r1, r2, r3, r4), (chi01, 1, 8, 1), require_integral=True
                    )
                except ValueError:
                    continue
                MM = explode(M, (alpha0, alpha1, alpha2, alpha3, alpha4))
                s = MM**-1 * MM.transpose()
                if (s - 1) ** 3 != 0:
                    continue
                data.append([[r0, r1, r2, r3, r4], [1, chi23, 1, chi40, chi01]])
    return data


def process_settings_stage_3():
    with open(f"{DATA_DIR}/5blocks_case_2_stage_2.txt") as f:
        lines = f.readlines()
        # To make uniform progress...
        lines = np.random.permutation(lines)
    with open(f"{DATA_DIR}/5blocks_case_2_stage_3.txt", "w") as f:
        for line in tqdm(lines, smoothing=0.0):
            s = ast.literal_eval(line)
            data = process_setting_stage_3(s)
            s["data"] = data
            if data != []:
                f.write(str(s) + "\n")
                f.flush()


def doall():
    # This takes 2000s on a 4 year old laptop
    # and ends with an empty file
    # 5blocks_case_2_stage_3.txt, meaning
    # that no solutions were found (as expected).
    #
    # Note: if you invoke this function, intermediate files,
    # generated during a previous run, will be overwritten.
    valid_settings()
    process_settings_stage_1()
    process_settings_stage_2()
    process_settings_stage_3()
    return


if __name__ == "__main__":
    t = time.time()
    doall()
    print(f"Test took {time.time() - t} seconds")
