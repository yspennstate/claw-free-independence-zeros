"""Checks for 'Two questions of Jerrum and Patel on zeros of the independence polynomial'.

Part 1  (Theorem 1.3) the graphs G_{m,h}: degree, transfer-matrix formula against a direct
        evaluation, absence of an induced S_{1,2,2} in G_{2,2}, and the zeros near 2^m e^{i pi/m}.
Part 2  (Theorem 1.1) a sanity check only: on random small claw-free graphs with random weights in
        R(2500(Delta+1)), Z_G != 0 and Re(Z_H/Z_B) > 0.97 at every vertex.
Part 3  (Proposition 1.2) the lower-bound examples: K_{Delta+1}, and the complements of C_5 and of
        the Wagner graph, with their independence polynomials computed from the graphs.

Requires Python 3 and mpmath.  Run:  python check.py
"""
import functools
import itertools
import math
import random

import mpmath as mp


# ---------------------------------------------------------------- general helpers
def indpoly(vertices, adj, weight):
    """Multivariate independence polynomial by deletion of the first vertex."""
    @functools.lru_cache(maxsize=None)
    def f(rem):
        if not rem:
            return 1
        v, rest = rem[0], rem[1:]
        return f(rest) + weight[v] * f(tuple(u for u in rest if u not in adj[v]))
    return f(tuple(sorted(vertices)))


def induced(adj, S):
    return {v: adj[v] & S for v in S}


def has_claw(adj):
    for v in adj:
        for a, b, c in itertools.combinations(adj[v], 3):
            if b not in adj[a] and c not in adj[a] and c not in adj[b]:
                return True
    return False


# ---------------------------------------------------------------- part 1
def build(m, h):
    """G_{m,h}: 2mh bags around a cycle; odd sites K2+K1, even sites 2K1, every m-th even site 3K1."""
    bags = []
    for _ in range(h):
        for q in range(m):
            bags.append("F")
            bags.append("S3" if q == m - 1 else "S2")
    adj, idx, n = {}, [], 0
    for t in bags:
        size = 2 if t == "S2" else 3
        idx.append(list(range(n, n + size)))
        n += size
    for v in range(n):
        adj[v] = set()
    for i, t in enumerate(bags):
        if t == "F":
            a, b, _ = idx[i]
            adj[a].add(b)
            adj[b].add(a)
        for x in idx[i]:
            for y in idx[(i + 1) % len(bags)]:
                adj[x].add(y)
                adj[y].add(x)
    return adj


def transfer(m, h, z):
    a, b, c = 2 * z * z + 3 * z, z * z + 2 * z, z ** 3 + 3 * z * z + 3 * z
    M = lambda x: mp.matrix([[1, x], [1, 0]])
    P = mp.eye(2)
    for _ in range(h):
        for q in range(m):
            P = P * M(a) * M(c if q == m - 1 else b)
    return P[0, 0] + P[1, 1]


def is_S122(adj, S):
    sub = induced(adj, S)
    deg = sorted(len(sub[v]) for v in S)
    if deg != [1, 1, 1, 2, 2, 3]:
        return False
    centre = next(v for v in S if len(sub[v]) == 3)
    seen, stack = {centre}, [centre]
    while stack:
        x = stack.pop()
        for y in sub[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    if len(seen) != 6:
        return False
    leaves = [v for v in S if len(sub[v]) == 1]
    return len({next(iter(sub[l])) for l in leaves}) == 3


def part1():
    mp.mp.dps = 30
    for m, h in [(2, 1), (3, 1), (2, 2), (4, 1)]:
        adj = build(m, h)
        assert len(adj) == h * (5 * m + 1)
        assert max(len(s) for s in adj.values()) <= 6
        for z in [mp.mpf(1), mp.mpc(0.3, -1.7), mp.mpf(-2), mp.mpc(-0.5, 0.25)]:
            direct = indpoly(adj.keys(), adj, {v: z for v in adj})
            tm = transfer(m, h, z)
            assert abs(direct - tm) < mp.mpf(10) ** -20 * max(1, abs(direct)), (m, h, z)
    print("part 1: orders h(5m+1), max degree <= 6, transfer formula = direct sum on (2,1),(3,1),(2,2),(4,1)")

    adj = build(2, 2)
    copies = sum(is_S122(adj, set(S)) for S in itertools.combinations(adj, 6))
    tree = {0: {1, 2, 4}, 1: {0}, 2: {0, 3}, 3: {2}, 4: {0, 5}, 5: {4}}
    assert is_S122(tree, set(tree))  # positive control: S_{1,2,2} itself
    assert copies == 0
    print("part 1: G_{2,2} (%d vertices): induced S_{1,2,2} found: %d (control on S_{1,2,2} itself: found)" % (len(adj), copies))

    mp.mp.dps = 80
    for m in [8, 12, 16, 20, 30]:
        def g(w):
            z = 2 ** m * w
            b = z * z + 2 * z
            return transfer(m, m, z) / (2 ** m * b ** m) ** m
        w = mp.findroot(g, mp.exp(1j * mp.pi / m))
        z = 2 ** m * w
        assert abs(g(w)) < mp.mpf(10) ** -40
        print("part 1: m=%2d  arg z_m = %s   pi/m = %s   |z_m|/2^m = %s   m^2|w_m - e^{i pi/m}| = %s"
              % (m, mp.nstr(mp.arg(z), 10), mp.nstr(mp.pi / m, 10), mp.nstr(abs(w), 8),
                 mp.nstr(m * m * abs(w - mp.exp(1j * mp.pi / m)), 4)))


# ---------------------------------------------------------------- part 2
def sample_R(K, rng):
    """A random point of R(K) = {(s+it)^2 : s >= 0, |t| < eta}, eta^2 = 1/(4K)."""
    eta = 1 / (2 * math.sqrt(K))
    s = rng.choice([rng.uniform(0, 10 * eta), rng.uniform(0, 3), rng.expovariate(0.05)])
    t = rng.uniform(-0.999, 0.999) * eta
    return complex(s, t) ** 2


def random_claw_free(n, p, rng):
    """The line graph of a random graph, or a random graph, kept if claw-free."""
    while True:
        if rng.random() < 0.5:
            k = max(3, n // 2)
            edges = [e for e in itertools.combinations(range(k), 2) if rng.random() < p]
            edges = edges[:n]
            adj = {i: set() for i in range(len(edges))}
            for i, j in itertools.combinations(range(len(edges)), 2):
                if set(edges[i]) & set(edges[j]):
                    adj[i].add(j)
                    adj[j].add(i)
        else:
            adj = {i: set() for i in range(n)}
            for i, j in itertools.combinations(range(n), 2):
                if rng.random() < p:
                    adj[i].add(j)
                    adj[j].add(i)
        if adj and not has_claw(adj):
            return adj


def part2(trials=300, seed=1):
    rng = random.Random(seed)
    worst = 10.0
    for _ in range(trials):
        adj = random_claw_free(rng.randint(4, 12), rng.uniform(0.2, 0.8), rng)
        Delta = max(len(s) for s in adj.values())
        K = 2500 * (Delta + 1)
        wt = {v: sample_R(K, rng) for v in adj}
        ZG = indpoly(adj.keys(), adj, wt)
        assert abs(ZG) > 0
        for v in adj:
            B = set(adj) - adj[v] - {v}
            ZB = indpoly(B, induced(adj, B), wt)
            ratio = ZG / ZB
            worst = min(worst, ratio.real)
    assert worst > 0.97
    print("part 2: %d random claw-free graphs, weights in R(2500(Delta+1)): min Re(Z_H/Z_B) over all vertices = %.4f > 0.97"
          % (trials, worst))


# ---------------------------------------------------------------- part 3
def coeffs(adj):
    """Coefficients of the univariate independence polynomial, by brute force."""
    n = len(adj)
    c = [0] * (n + 1)
    for mask in range(1 << n):
        S = [v for v in range(n) if mask >> v & 1]
        if all(b not in adj[a] for a, b in itertools.combinations(S, 2)):
            c[len(S)] += 1
    while c and c[-1] == 0:
        c.pop()
    return c


def part3():
    for Delta in range(1, 8):
        n = Delta + 1
        adj = {v: set(range(n)) - {v} for v in range(n)}
        c = coeffs(adj)
        assert c == [1, n] and not has_claw(adj)
        assert max(len(s) for s in adj.values()) == Delta
        root = -c[0] / c[1]
        assert abs(root + 1 / (Delta + 1)) < 1e-12
    # complements of triangle-free graphs on N = r(3,k+1)-1 vertices with no independent set of size k+1:
    # C_5 (N = 5, k = 2) and the Wagner graph, the circulant C_8(1,4) (N = 8, k = 3)
    tests = [(5, [1]), (8, [1, 4])]
    for N, steps in tests:
        T = {v: {(v + s) % N for s in steps} | {(v - s) % N for s in steps} for v in range(N)}
        assert all(not (T[a] & T[b]) for a in T for b in T[a])  # triangle-free
        G = {v: set(range(N)) - T[v] - {v} for v in range(N)}  # complement
        assert not has_claw(G)
        c = coeffs(G)
        assert len(c) == 3 and c[1] == N and c[2] == len([1 for a in T for b in T[a] if a < b])
        a2 = c[2]
        assert 4 * a2 <= N * N
        root = max(r.real for r in mp.polyroots(c[::-1]))
        assert -2.0 / N <= root <= -1.0 / N, (N, root)
        print("part 3: complement of a triangle-free graph on %d vertices: Z = %s, zero nearest 0 at %.4f in [-2/N,-1/N]"
              % (N, c, root))
    print("part 3: K_{Delta+1}, Delta=1..7: Z = 1+(Delta+1)z, zero at -1/(Delta+1), claw-free, max degree Delta")


if __name__ == "__main__":
    part1()
    part2()
    part3()
    print("all checks passed")
