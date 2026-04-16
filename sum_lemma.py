from functools import cache


@cache
def sum_lemma_max(g, a):
    """
    This mplements the algorithm from the proof of Lemma 13.3.
    Assume g_1(x),..,g_n(x) x>=1 are descending with limit 0.
    For a given a we compute the maximum value
    for g_1(x_1)+...+g_n(x_n) under the condition that
    this sum is < a.
    """
    assert a > 0
    n = len(g)
    f = g[0]
    if n == 1:
        x = 1
        while f(x) >= a:
            x += 1
        return f(x)
    else:
        m = sum_lemma_max(g[1:], a)

        def m_(x):
            return sum_lemma_max(g[1:], a - f(x))

        x = 1
        first_term = -1
        while f(x) >= a - m:
            if f(x) >= a:
                x += 1
                continue
            if f(x) + m_(x) > first_term:
                first_term = f(x) + m_(x)
            x += 1
        second_term = sum_lemma_max((f,), a - m) + m
        return max(first_term, second_term)


@cache
def sum_lemma_total(g, a):
    """
    This implements the algorithm from the proof of Lemma 13.3.
    Assume g_1(x),..,g_n(x) x>=1 are descending with limit 0.
    For a given a we compute the solutions of
    g_1(x_1)+...+g_n(x_n) = a.
    """
    assert a > 0
    if len(g) == 0:
        return ()
    f = g[0]

    def finv(a):
        x = 1
        while f(x) > a:
            x += 1
        if f(x) == a:
            return x
        else:
            return None

    if len(g) == 1:
        if (x := finv(a)) is not None:
            return ((x,),)
        else:
            return ()
    m = sum_lemma_max(g[1:], a)
    x = 1
    sols = []
    while f(x) >= a - m:
        if f(x) >= a:
            x += 1
            continue
        sols_new = sum_lemma_total(g[1:], a - f(x))
        sols.extend((x,) + s for s in sols_new)
        x += 1
    return tuple(sols)


if __name__ == "__main__":
    from fractions import Fraction
    from itertools import product

    from tqdm import tqdm

    def f(n):
        return Fraction(1, n)

    sols = sum_lemma_total((f, f, f, f), 1)
    naive_sols = []
    for p in tqdm(product(range(1, 50), repeat=4), total=49**4):  # 50 is a guess
        if sum([f(p[i]) for i in range(0, 4)]) == 1:
            naive_sols.append(p)
    assert set(sols) == set(naive_sols)
