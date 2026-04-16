import pprint
from fractions import Fraction

from sage.all import ZZ, Integer, sqrt

from mutations import Quiver, explode, gram_matrix, implode
from sum_lemma import sum_lemma_total


def valid_settings():
    sols = []
    for chi12 in range(1, 3):
        for alpha1 in range(1, 5):
            for alpha2 in range(1, 5):
                if chi12**2 * alpha1 * alpha2 != 4:
                    continue
                for chi30 in range(1, 3):
                    for alpha0 in range(1, 5):
                        for alpha3 in range(1, 5):
                            if chi30**2 * alpha3 * alpha0 != 4:
                                continue
                            s = alpha0 + alpha1 + alpha2 + alpha3
                            if s > 11:
                                continue
                            K2 = 12 - s
                            sols.append(
                                (K2, chi12, chi30, alpha0, alpha1, alpha2, alpha3)
                            )
    sols = list(sorted(sols))
    return sols


def find_gram_matrices(sol):
    sols = []

    K2 = sol[0]
    chi12 = sol[1]
    chi30 = sol[2]
    alpha0 = sol[3]
    alpha1 = sol[4]
    alpha2 = sol[5]
    alpha3 = sol[6]

    def f0(r0):
        return Fraction(alpha1, r0 * r0)

    def f1(r1):
        return Fraction(alpha0, r1 * r1)

    a = Fraction(K2 * alpha0 * alpha1, 4)

    rank_pairs = sum_lemma_total((f0, f1), a)

    for r0, r1 in rank_pairs:
        r2 = sqrt(Integer(alpha1) * r1**2 / alpha2)
        if not r2.is_integer():
            continue
        r2 = ZZ(r2)
        r3 = sqrt(Integer(alpha0) * r0**2 / alpha3)
        if not r3.is_integer():
            continue
        r3 = ZZ(r3)

        gap = K2 - Integer(chi12) / (r1 * r2) - Integer(chi30) / (r3 * r0)
        # We have gap == chi23 / (r2 * r3)+ chi01 / (r0 * r1)
        chi01_upperbound = (gap * r0 * r1).floor()
        chi23_upperbound = (gap * r2 * r3).floor()

        for chi01 in range(1, chi01_upperbound + 1):
            for chi23 in range(1, chi23_upperbound + 1):
                if gap != Integer(chi23) / (r2 * r3) + Integer(chi01) / (r0 * r1):
                    continue

                # check that the origin is in the forbidden area
                # some of these conditions are no doubt redundant
                if 2 * r0 > chi30 * alpha3 * r3:
                    continue
                if 2 * r3 > chi30 * alpha0 * r0:
                    continue
                if 2 * r1 > chi01 * alpha0 * r0:
                    continue
                if 2 * r2 > chi23 * alpha3 * r3:
                    continue
                try:
                    G = gram_matrix(
                        [r0, r1, r2, r3], [chi01, chi12, chi23], require_integral=True
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

                # check that the collection is block complete

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
                        "alphas": [
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
