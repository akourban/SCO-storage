# naive.py
import sys

# La fonction size_2D de dyprosco pour définir la taille de notre tableau donc 
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
L= [] # Listes de toutes les arêtes internes(atomes voisins sur la grille n1xn2)


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

sizeC = max(n2 + 2 * n1, len(B) + 1)  # Borne supérieure pour l'indice du bord c
smax = len(L) # Nombre maximales de liaisons alignées sur un multiple de 4 
if smax % 4:
    smax += 4 - smax % 4
size2D = size_2D(n1, n2, sizeC)

d = [0] * size2D # Tableau 1D de densité d'états

for i in range(2**n):
    bits = [(i >> b) & 1 for b in range(n)]               # état de la molécule, 0 ou 1

    m = sum(bits)                                         # atomes à 1
    s = sum(1 for a, b in L if bits[a] != bits[b])        # arêtes déséquilibrées
    c = sum(bits[a] for a in B)                           # atomes à 1 sur le bord

    k = index_2D(m, s, c, n, sizeC, smax, size2D)
    if k < 0 or k >= size2D:
        sys.exit(f"index hors tableau : i={i} (m={m}, s={s}, c={c}) -> {k}")
    d[k] += 1

for m in range(n + 1):
    for s in list(range(smax+1)):
        for c in range(sizeC):
            k = index_2D(m, s, c, n, sizeC, smax, size2D)
            if 0 <= k < size2D and d[k] > 0:
                print(f"{m} {s} {c} {d[k]}")

## Afin de comparer directement les résultats bruts avec le code C de dyprosco,
## nous avons adapté write_2dCart_msc() dans src/d-io.c pour exporter les indices
## entiers positifs (m, s, c) et la dégénérescence associée, plutôt que d'appliquer
## la conversion linéaire vers les grandeurs physiques signées (M, E_interne, E_bord).


print("\nVisualisation de l'espace alloué (triangle/losange) :")
header = " m \\ s |" + "".join(f"{s:>5}" for s in range(smax + 1))
print(header)
print("-" * len(header))

for m in range(n + 1):
    row = f"{m:>5}  |"
    for s in range(smax + 1):
        # On teste c=0 pour savoir si le couple (m, s) est alloué dans le stockage creux
        k0 = index_2D(m, s, 0, n, sizeC, smax, size2D)
        if k0 < 0 or k0 >= size2D:
            row += f"{'.':>5}"  # Hors du losange (non alloué en mémoire)
        else:
            # Somme sur toute la tranche c pour ce (m, s) alloué
            total_ms = sum(d[k0 + c] for c in range(sizeC))
            row += f"{total_ms:>5}"
    print(row)