import pprint
from fractions import Fraction

from sage.all import Integer, Rational, ceil

from mutations import Quiver, explode, gram_matrix, implode


def valid_settings():
    sols = []
    chi12 = 1
    for alpha1 in range(1, 4):
        for alpha2 in range(1, 4):
            p1 = chi12**2 * alpha1 * alpha2
            if p1 > 3:  # (13.38)
                continue
            for chi30 in range(1, 5):  # (13.30)
                # max value of alpha_i = 11-3=8
                for alpha0 in range(1, 9):
                    for alpha3 in range(1, 9):
                        s = alpha0 + alpha1 + alpha2 + alpha3
                        if s > 11:
                            continue
                        K2 = 12 - s
                        p2 = chi30**2 * alpha3 * alpha0
                        if p2 < 5:
                            continue
                        if p1 * p2 > 16:
                            continue
                        sols.append((K2, chi12, chi30, alpha0, alpha1, alpha2, alpha3))
    sols = list(sorted(sols))
    return sols


def find_gram_matrices(sol):
    K2 = sol[0]
    chi12 = sol[1]
    chi30 = sol[2]
    alpha0 = sol[3]
    alpha1 = sol[4]
    alpha2 = sol[5]
    alpha3 = sol[6]
    sols = []
    r03_upperbound = 22 * chi30  # (13.46)
    for r0 in range(1, r03_upperbound + 1):
        for r3 in range(1, r03_upperbound + 1):
            # forbidden area check (13.34)
            if 2 * r0 > r3 * chi30 * alpha3:
                continue
            if 2 * r3 > r0 * chi30 * alpha0:
                continue

            gap = (
                K2
                - Fraction(2, (alpha0 * r0**2))
                - Fraction(2, (alpha3 * r3**2))
                - Fraction(chi30, (r0 * r3))
            )  # (13.34)
            if gap <= 0:
                continue

            r1r2_lower = ceil(1 / Rational(gap))

            gamma = 1 / (Fraction(1, 2) - Fraction(2, (chi30**2 * alpha3 * alpha0)))
            # (13.44)
            if (
                1 / r1r2_lower
                + Rational(gamma) / (alpha0 * r0**2)
                + Rational(gamma) / (alpha3 * r3**2)
                + Integer(chi30) / (r0 * r3)
                < K2
            ):
                continue

            sq = alpha0 * r0**2 + alpha3 * r3**2
            for r1 in range(1, sq + 1):
                if alpha1 * r1**2 > sq:
                    continue
                for r2 in range(1, sq + 1):
                    if alpha1 * r1**2 + alpha2 * r2**2 > sq:
                        continue
                    if r1 * r2 < r1r2_lower:
                        continue

                    gap2 = K2 - Integer(chi12) / (r1 * r2) - Integer(chi30) / (r3 * r0)
                    # We have gap2 == chi23 / (r2 * r3)+ chi01 / (r0 * r1)
                    chi01_upperbound = (gap2 * r0 * r1).floor()
                    chi23_upperbound = (gap2 * r2 * r3).floor()

                    for chi01 in range(1, chi01_upperbound + 1):
                        for chi23 in range(1, chi23_upperbound + 1):
                            if gap2 != Integer(chi23) / (r2 * r3) + Integer(chi01) / (
                                r0 * r1
                            ):
                                continue
                            # forbidden area check
                            if 2 * r1 > r0 * chi01 * alpha0:
                                continue
                            if 2 * r2 > r3 * chi23 * alpha3:
                                continue
                            try:
                                G = gram_matrix(
                                    [r0, r1, r2, r3],
                                    [chi01, chi12, chi23],
                                    require_integral=True,
                                )
                            except ValueError:
                                continue

                            alphas = [alpha0, alpha1, alpha2, alpha3]
                            GE = explode(G, alphas)
                            GEinv = GE**-1
                            if (GE - GE.transpose()).rank() > 2:
                                continue
                            if (GEinv * GE.transpose() - 1) ** 3 != 0:
                                continue

                            Q = Quiver(GE)
                            if not Q.is_block_complete():
                                continue

                            # make sure that the vertex (1,2) is not-admissible
                            Ginv = implode(GEinv, alphas)
                            if Ginv[0, 2] > 0:
                                continue
                            if Ginv[1, 3] > 0:
                                continue

                            sols.append(
                                {
                                    "K2": K2,
                                    "r": [r0, r1, r2, r3],
                                    "alpha": [
                                        alpha0,
                                        alpha1,
                                        alpha2,
                                        alpha3,
                                    ],
                                    "Mred": [list(l) for l in G.rows()],
                                }
                            )
    return sols


def full_list():
    gms = []
    sols = valid_settings()
    for sol in sols:
        gm = find_gram_matrices(sol)
        gms.extend(gm)
    return gms


if __name__ == "__main__":
    gms = full_list()
    print(f"{len(gms)} solutions found!\n")
    pprint.pprint(gms)
