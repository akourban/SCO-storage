# naive.py
import sys


def size_2D(n1, n2, sizeC):
    n = n1 * n2
    if sizeC < n2 + 2 * n1:
        sizeC = n2 + 2 * n1
    smax = n1 * (n2 - 1) + n2 * (n1 - 1)
    rs = smax % 4
    if rs != 0:
        smax += 4 - rs
    if n % 2 == 0:
        ll = n // 2 - 1
        return 4 * ll * (ll + 1) * sizeC + smax * sizeC + (n + 1) * sizeC
    else:
        ll = (n - 1) // 2 - 1
        return 4 * ll * (ll + 1) * sizeC + 2 * smax * sizeC + (n + 1) * sizeC


def index_2D(m, s, c, n, sizeC, smax, size2D):
    if s == 0:
        return size2D - (n + 1) * sizeC + m * sizeC + c

    if n % 2 == 0:
        ll = n // 2 - 1
        sizeP1 = 2 * ll * (ll + 1) * sizeC
        if m < n // 2 + 1:
            if s <= 4 * m and s <= smax:
                return 2 * (m - 1) * m * sizeC + (s - 1) * sizeC + c
            return -1
        ii = m - 1 - n // 2
        kk = 2 * (ll * (ll + 1) - (ll - ii) * (ll - ii + 1))
        if s <= 4 * (n - m):
            return sizeP1 + smax * sizeC + kk * sizeC + (s - 1) * sizeC + c
        return -1
    else:
        ll = (n - 1) // 2 - 1
        sizeP1 = 2 * ll * (ll + 1) * sizeC + smax * sizeC
        if m <= (n - 1) // 2:
            if s <= 4 * m and s <= smax:
                return 2 * (m - 1) * m * sizeC + (s - 1) * sizeC + c
            return -1
        elif m == (n - 1) // 2 + 1:
            if s <= smax:
                return sizeP1 + (s - 1) * sizeC + c
            return -1
        else:
            ii = m - 2 - (n - 1) // 2
            kk = 2 * (ll * (ll + 1) - (ll - ii) * (ll - ii + 1))
            if s <= 4 * (n - m):
                return sizeP1 + smax * sizeC + kk * sizeC + (s - 1) * sizeC + c
            return -1



n1, n2 = (3, 3)
n = n1 * n2

L = []

for r in range(n1):
    for col in range(n2):
        k = r * n2 + col
        if col < n2 - 1:
            L.append((k, k + 1))
        if r < n1 - 1:
            L.append((k, k + n2))

# B : atomes du bord
B = [r * n2 + col
     for r in range(n1) for col in range(n2)
     if r in (0, n1 - 1) or col in (0, n2 - 1)]

# paramètres de stockage (mêmes que size_2D)
sizeC = max(n2 + 2 * n1, len(B) + 1)     # il faut c dans [0, |B|]
smax = len(L)
if smax % 4:
    smax += 4 - smax % 4
size2D = size_2D(n1, n2, sizeC)

d = [0] * size2D

for i in range(2**n):
    bits = [(i >> b) & 1 for b in range(n)]               # la molécule

    m = sum(bits)                                         # atomes à 1
    s = sum(1 for a, b in L if bits[a] != bits[b])        # arêtes déséquilibrées
    c = sum(bits[a] for a in B)                           # atomes à 1 sur le bord

    k = index_2D(m, s, c, n, sizeC, smax, size2D)
    if k < 0 or k >= size2D:
        sys.exit(f"index hors tableau : i={i} (m={m}, s={s}, c={c}) -> {k}")
    d[k] += 1

print(f"size2D = {size2D}, total = {sum(d)} (attendu {2**n})")

# ---------- affichage de quelques molécules avec leurs m, s, c ----------

def calcul_msc(i):
    bits = [(i >> b) & 1 for b in range(n)]
    m = sum(bits)
    s = sum(1 for a, b in L if bits[a] != bits[b])
    c = sum(bits[a] for a in B)
    return bits, m, s, c

def afficher(indices, largeur=16):
    blocs = []
    for i in indices:
        bits, m, s, c = calcul_msc(i)
        bloc = [f"i = {i}"]
        for r in range(n1):
            bloc.append(" ".join(str(bits[r * n2 + col]) for col in range(n2)))
        bloc.append(f"m = {m}")
        bloc.append(f"s = {s}")
        bloc.append(f"c = {c}")
        blocs.append(bloc)

    for lignes in zip(*blocs):
        print("".join(l.ljust(largeur) for l in lignes))


def afficher_triangle():
    print("  m \\ s " + "".join(f"{s:6d}" for s in range(smax + 1)))
    for m in range(n + 1):
        ligne = f"  {m:2d}    "
        for s in range(smax + 1):
            if index_2D(m, s, 0, n, sizeC, smax, size2D) < 0:
                ligne += "     ."
            else:
                tot = sum(d[index_2D(m, s, c, n, sizeC, smax, size2D)]
                          for c in range(len(B) + 1))
                ligne += f"{tot:6d}"
        print(ligne)

afficher_triangle()