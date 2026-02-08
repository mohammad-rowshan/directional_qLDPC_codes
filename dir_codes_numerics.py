# Copyright (c) 2026 Mohammad Rowshan
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to use,
# copy, modify, and distribute the Software for **scientific and academic
# research purposes only**, subject to the following conditions:
#
# 1. Commercial use of this Software, including but not limited to use in
#    proprietary products, services, or for-profit activities, is strictly
#    prohibited without prior written permission from the copyright holder.
#
# 2. The above copyright notice and this permission notice shall be included in
#    all copies or substantial portions of the Software.
#
# 3. The Software is provided "as is", without warranty of any kind, express or
#    implied, including but not limited to the warranties of merchantability,
#    fitness for a particular purpose, and noninfringement.
#
# This notice is a permissive research-only license inspired by the MIT License.


"""
Numerics for "Directional (route-generated) qLDPC codes" in the square-grid checkerboard model.

What it does (static code properties):
- Parse direction words (e.g. "N2E3N2", "NE^2NE^2N", "NE2NE2N")
- Compute route offsets Q_j = S_{j-1} + S_j (Lemma route-to-offset in the manuscript)
- Build CSS parity-check matrices H_X, H_Z on a checkerboard torus Z_{Lx} x Z_{Ly}
- Validate commutation H_X H_Z^T = 0 (overlap parity even)
- Compute (n,k) exactly via GF(2) rank
- Compute small-weight distances exactly up to a chosen cutoff by exhaustive search
- Optionally enumerate coset-constant layouts derived from Δ_odd(W)
- Optionally scan over words (exhaustive) for a fixed torus and report best candidates

Design goal:
- Run by pressing F5 in Spyder (no CLI / argparse needed).
- Print results as plain tables in the console.
- Optional CSV export (OFF by default).
- No LaTeX table generation.

Notes:
- Distance search is combinatorial; keep DIST_MAX small (e.g., 3–5 for n=64; 2–3 for n>150).
- This script focuses on the square-grid (N,E,S,W) instantiation used in the manuscript’s main case study.
"""

from __future__ import annotations

import csv
import itertools
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# =============================================================================
# 0) Helpers
# =============================================================================

def popcount(x: int) -> int:
    """Return number of 1-bits (Python 3.8+)."""
    try:
        return x.bit_count()  # type: ignore[attr-defined]
    except AttributeError:
        return bin(x).count("1")


def print_table(rows: List[Dict[str, object]], columns: List[str], title: str = "", max_rows: int = 50) -> None:
    """
    Print a nice fixed-width table to the console.

    rows: list of dicts (same keys)
    columns: column order to print
    """
    if title:
        print("\n" + title)
        print("-" * len(title))

    if not rows:
        print("(no rows)")
        return

    # limit rows
    shown = rows[:max_rows]

    # compute column widths
    def cell_str(val: object) -> str:
        return "" if val is None else str(val)

    widths: Dict[str, int] = {}
    for c in columns:
        widths[c] = max(len(c), *(len(cell_str(r.get(c, ""))) for r in shown))

    # header
    header = "  ".join(c.ljust(widths[c]) for c in columns)
    print(header)
    print("  ".join("-" * widths[c] for c in columns))

    # body
    for r in shown:
        line = "  ".join(cell_str(r.get(c, "")).ljust(widths[c]) for c in columns)
        print(line)

    if len(rows) > max_rows:
        print(f"... ({len(rows) - max_rows} more rows not shown)")


def maybe_write_csv(path: str, rows: List[Dict[str, object]]) -> None:
    """Optional CSV export."""
    if not rows:
        return
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"[saved] {path} ({len(rows)} rows)")


# =============================================================================
# 1) Parsing / geometry for square-grid direction words
# =============================================================================

DIR_VECS: Dict[str, Tuple[int, int]] = {
    "N": (0, 1),
    "E": (1, 0),
    "S": (0, -1),
    "W": (-1, 0),
}

OPPOSITE: Dict[str, str] = {"N": "S", "S": "N", "E": "W", "W": "E"}


def parse_direction_word(s: str) -> List[str]:
    """
    Parse words like 'N2E3N2', 'NE^2NE^2N', 'NE2NE2N' into an expanded list of direction letters.
    Non-direction characters are ignored (useful if pasting LaTeX snippets).
    """
    s = (s or "").strip().replace(" ", "")
    out: List[str] = []
    i = 0
    while i < len(s):
        ch = s[i]
        if ch in "NESW":
            letter = ch
            i += 1
            rep: Optional[int] = None

            # exponent forms: ^2 or ^{2}
            if i < len(s) and s[i] == "^":
                i += 1
                if i < len(s) and s[i] == "{":
                    i += 1
                    j = i
                    while j < len(s) and s[j].isdigit():
                        j += 1
                    if j == i:
                        raise ValueError(f"Expected digits after '^{{' in '{s}'")
                    rep = int(s[i:j])
                    if j >= len(s) or s[j] != "}":
                        raise ValueError(f"Expected '}}' after exponent in '{s}'")
                    i = j + 1
                else:
                    j = i
                    while j < len(s) and s[j].isdigit():
                        j += 1
                    if j == i:
                        raise ValueError(f"Expected digits after '^' in '{s}'")
                    rep = int(s[i:j])
                    i = j
            else:
                # direct repetition: E2, N10, ...
                j = i
                while j < len(s) and s[j].isdigit():
                    j += 1
                if j > i:
                    rep = int(s[i:j])
                    i = j

            if rep is None:
                rep = 1
            out.extend([letter] * rep)
        else:
            i += 1

    if not out:
        raise ValueError(f"No directions parsed from '{s}'")
    return out


def compress_word(word: Sequence[str]) -> str:
    """Run-length encode (lossless) to a compact 'N2E3N2'-style string."""
    if not word:
        return ""
    parts: List[str] = []
    cur = word[0]
    cnt = 1
    for ch in word[1:]:
        if ch == cur:
            cnt += 1
        else:
            parts.append(cur + (str(cnt) if cnt != 1 else ""))
            cur, cnt = ch, 1
    parts.append(cur + (str(cnt) if cnt != 1 else ""))
    return "".join(parts)


def route_offsets(word_dirs: Sequence[str]) -> List[Tuple[int, int]]:
    """
    Compute offsets Q_j = S_{j-1}+S_j for the expanded word (square-grid).
    Offsets are integer vectors in Z^2.
    """
    S_prev = (0, 0)
    S = (0, 0)
    offsets: List[Tuple[int, int]] = []
    for d in word_dirs:
        dx, dy = DIR_VECS[d]
        S = (S[0] + dx, S[1] + dy)
        Q = (S_prev[0] + S[0], S_prev[1] + S[1])
        offsets.append(Q)
        S_prev = S
    return offsets


def support_weight_mod2(offsets: Sequence[Tuple[int, int]]) -> int:
    """
    Stabilizer support is over F2: repeated offsets cancel.
    Returns the number of offsets with odd multiplicity.
    """
    counts: Dict[Tuple[int, int], int] = {}
    for q in offsets:
        counts[q] = counts.get(q, 0) ^ 1
    return sum(counts.values())


def no_immediate_backtrack(word_dirs: Sequence[str]) -> bool:
    """Reject patterns like '...N S...' or '...E W...'."""
    for a, b in zip(word_dirs, word_dirs[1:]):
        if OPPOSITE[a] == b:
            return False
    return True


def offsets_distinct(offsets: Sequence[Tuple[int, int]]) -> bool:
    """Require all offsets Q_j distinct (no self-hit)."""
    return len(set(offsets)) == len(offsets)


# =============================================================================
# 2) Odd-multiplicity differences Δ_odd(W) (layout invariant)
# =============================================================================

def delta_odd_set(offsets: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Δ_odd(W) = {v : multiplicity of v among pairwise differences Q_j-Q_i (i<j) is odd}.
    Returns a (not-necessarily-minimal) generator set; v=(0,0) omitted.
    """
    mu: Dict[Tuple[int, int], int] = {}
    m = len(offsets)
    for i in range(m):
        xi, yi = offsets[i]
        for j in range(i + 1, m):
            xj, yj = offsets[j]
            v = (xj - xi, yj - yi)
            mu[v] = mu.get(v, 0) ^ 1
    return [v for (v, parity) in mu.items() if parity == 1 and v != (0, 0)]


# =============================================================================
# 3) Torus geometry + layouts + parity-check matrices
# =============================================================================

def coords_data_ancillas(Lx: int, Ly: int) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
    data: List[Tuple[int, int]] = []
    anc: List[Tuple[int, int]] = []
    for y in range(Ly):
        for x in range(Lx):
            if (x + y) % 2 == 0:
                data.append((x, y))
            else:
                anc.append((x, y))
    return data, anc


def layout_row_alternating(_Lx: int, _Ly: int):
    """Example layout used in the manuscript: X on even y, Z on odd y (on ancilla sites)."""
    def alpha(a: Tuple[int, int]) -> int:
        _x, y = a
        return 0 if (y % 2 == 0) else 1  # 0=X, 1=Z
    return alpha


class UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, a: int) -> int:
        while self.parent[a] != a:
            self.parent[a] = self.parent[self.parent[a]]
            a = self.parent[a]
        return a

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def cosets_from_generators(
    Lx: int, Ly: int, ancillas: Sequence[Tuple[int, int]], gens: Sequence[Tuple[int, int]]
) -> Tuple[Dict[Tuple[int, int], int], int]:
    """
    Compute cosets of the subgroup generated by 'gens' acting on ancilla coordinates modulo (Lx, Ly).
    Returns (coset_id_map, num_cosets).
    """
    idx = {a: i for i, a in enumerate(ancillas)}
    uf = UnionFind(len(ancillas))

    for a in ancillas:
        ia = idx[a]
        for gx, gy in gens:
            b = ((a[0] + gx) % Lx, (a[1] + gy) % Ly)
            if b in idx:
                uf.union(ia, idx[b])
            b = ((a[0] - gx) % Lx, (a[1] - gy) % Ly)
            if b in idx:
                uf.union(ia, idx[b])

    root_to_cid: Dict[int, int] = {}
    coset_id: Dict[Tuple[int, int], int] = {}
    for a in ancillas:
        r = uf.find(idx[a])
        if r not in root_to_cid:
            root_to_cid[r] = len(root_to_cid)
        coset_id[a] = root_to_cid[r]
    return coset_id, len(root_to_cid)


def layout_from_cosets(ancillas: Sequence[Tuple[int, int]], coset_id: Dict[Tuple[int, int], int], assignment: Sequence[int]):
    """Turn an assignment on cosets (0=X,1=Z) into a layout function alpha(a)."""
    def alpha(a: Tuple[int, int]) -> int:
        return int(assignment[coset_id[a]])
    return alpha


def build_H(Lx: int, Ly: int, word_dirs: Sequence[str], alpha) -> Tuple[List[int], List[int], int]:
    """
    Build H_X and H_Z as lists of bit-rows (Python int), where bit i refers to the i-th data qubit.
    Returns (Hx_rows, Hz_rows, n_data).
    """
    if Lx % 2 or Ly % 2:
        raise ValueError("Require Lx, Ly even for checkerboard partition.")
    data, anc = coords_data_ancillas(Lx, Ly)
    n = len(data)
    data_index = {q: i for i, q in enumerate(data)}
    offsets = route_offsets(word_dirs)

    Hx: List[int] = []
    Hz: List[int] = []
    for a in anc:
        row = 0
        for Qx, Qy in offsets:
            q = ((a[0] + Qx) % Lx, (a[1] + Qy) % Ly)
            if q not in data_index:
                # If this triggers, your coordinate convention is inconsistent for this (Lx,Ly,word).
                raise ValueError(f"Parity mismatch: anchor={a}, offset={(Qx,Qy)}, target={q}")
            row ^= (1 << data_index[q])  # F2 support (duplicates cancel)
        if alpha(a) == 0:
            Hx.append(row)
        else:
            Hz.append(row)
    return Hx, Hz, n


# =============================================================================
# 4) GF(2) linear algebra on bit-rows
# =============================================================================

def gf2_rank(rows: Sequence[int]) -> Tuple[int, Dict[int, int]]:
    """
    Return (rank, basis) where basis maps pivot_bit_position -> row.
    """
    basis: Dict[int, int] = {}
    rank = 0
    for r in rows:
        x = int(r)
        while x:
            p = x.bit_length() - 1
            if p in basis:
                x ^= basis[p]
            else:
                basis[p] = x
                rank += 1
                break
    return rank, basis


def gf2_in_rowspace(v: int, basis: Dict[int, int]) -> bool:
    """Check if v is in rowspace spanned by basis."""
    x = int(v)
    while x:
        p = x.bit_length() - 1
        if p in basis:
            x ^= basis[p]
        else:
            return False
    return True


def commute_ok(Hx: Sequence[int], Hz: Sequence[int]) -> bool:
    """Check CSS commutation: all overlaps have even parity."""
    for rx in Hx:
        for rz in Hz:
            if popcount(rx & rz) & 1:
                return False
    return True


def build_column_syndromes(rows: Sequence[int], n: int) -> List[int]:
    """
    For check rows (m rows), return column syndrome masks (m-bit ints):
    col_syn[i] has bit j set iff row j includes column i.
    Then syndrome(subset) = XOR of col_syn over i in subset.
    """
    col = [0] * n
    for j, r in enumerate(rows):
        x = int(r)
        while x:
            lsb = x & -x
            i = lsb.bit_length() - 1
            col[i] |= (1 << j)
            x ^= lsb
    return col


def distance_search_upto(H_comm: Sequence[int], other_basis: Dict[int, int], n: int, w_max: int) -> Optional[int]:
    """
    Exhaustive search for the minimum weight <= w_max of a nontrivial logical operator:
      - kernel condition: commute with all checks in H_comm (syndrome==0)
      - not in rowspace of other type (other_basis)
    Returns exact distance if <= w_max; else None meaning d > w_max.
    """
    col_syn = build_column_syndromes(H_comm, n)
    for w in range(1, w_max + 1):
        for comb in itertools.combinations(range(n), w):
            syn = 0
            v = 0
            for i in comb:
                syn ^= col_syn[i]
                v |= 1 << i
            if syn != 0:
                continue
            if gf2_in_rowspace(v, other_basis):
                continue
            return w
    return None


# =============================================================================
# 5) Evaluation results
# =============================================================================

@dataclass
class CodeResult:
    word: str
    word_len: int
    support_wt: int
    Lx: int
    Ly: int
    layout: str
    commute: bool
    n: int
    k: int
    rankX: int
    rankZ: int
    dX: Optional[int]
    dZ: Optional[int]
    num_cosets: Optional[int]


def fmt_d(d: Optional[int], dist_max: int) -> str:
    """Pretty distance cell for console."""
    if dist_max <= 0:
        return "-"
    if d is None:
        return f">{dist_max}"
    return str(d)


def evaluate_word(
    word_str: str,
    Lx: int,
    Ly: int,
    layout_mode: str = "row",        # "row", "coset", or "row+coset"
    dist_max: int = 0,
    max_coset_layouts: int = 256,    # safety cap for coset enumeration
) -> List[CodeResult]:
    """
    Evaluate a word on a given torus for one or more layouts.
    Returns a list of CodeResult, one per layout tested.
    """
    word_dirs = parse_direction_word(word_str)
    offsets = route_offsets(word_dirs)
    support_wt = support_weight_mod2(offsets)

    # generators for coset layouts
    gens = delta_odd_set(offsets)

    results: List[CodeResult] = []
    tested_layouts: List[Tuple[str, object, Optional[int]]] = []

    # Always include row layout if requested
    if layout_mode in ("row", "row+coset"):
        alpha = layout_row_alternating(Lx, Ly)
        tested_layouts.append(("row", alpha, None))

    # Coset layouts (constant on cosets of L(W))
    if layout_mode in ("coset", "row+coset"):
        _data, anc = coords_data_ancillas(Lx, Ly)
        coset_id, c = cosets_from_generators(Lx, Ly, anc, gens)

        # enumerate all assignments up to global flip: fix coset 0 to X (0)
        num_assignments = 2 ** max(c - 1, 0)
        if num_assignments > max_coset_layouts:
            # print a warning and skip (or you can raise)
            print(f"[warn] coset enumeration too large: 2^{c-1} = {num_assignments} > {max_coset_layouts}. Skipping coset layouts.")
        else:
            for bits in itertools.product([0, 1], repeat=max(c - 1, 0)):
                assignment = [0] + list(bits)
                alpha = layout_from_cosets(anc, coset_id, assignment)
                tested_layouts.append((f"coset{assignment}", alpha, c))

    # evaluate each tested layout
    for layout_name, alpha, num_cosets in tested_layouts:
        try:
            Hx, Hz, n = build_H(Lx, Ly, word_dirs, alpha)
        except ValueError:
            continue

        comm = commute_ok(Hx, Hz)
        if not comm:
            continue

        rx, bx = gf2_rank(Hx)
        rz, bz = gf2_rank(Hz)
        k = n - rx - rz

        dX = dZ = None
        if dist_max > 0 and k > 0:
            dX = distance_search_upto(Hz, bx, n, dist_max)  # X-logicals commute with Z checks
            dZ = distance_search_upto(Hx, bz, n, dist_max)  # Z-logicals commute with X checks

        results.append(
            CodeResult(
                word=compress_word(word_dirs),
                word_len=len(word_dirs),
                support_wt=support_wt,
                Lx=Lx,
                Ly=Ly,
                layout=layout_name,
                commute=comm,
                n=n,
                k=k,
                rankX=rx,
                rankZ=rz,
                dX=dX,
                dZ=dZ,
                num_cosets=num_cosets,
            )
        )

    return results


# =============================================================================
# 6) Symmetry reduction (optional, for word searches)
# =============================================================================

def apply_dir_map(word: Sequence[str], mapping: Dict[str, str]) -> List[str]:
    return [mapping[ch] for ch in word]


def reversed_inverted(word: Sequence[str]) -> List[str]:
    return [OPPOSITE[ch] for ch in reversed(word)]


D4_MAPS: List[Dict[str, str]] = [
    # rotations
    {"N": "N", "E": "E", "S": "S", "W": "W"},  # id
    {"N": "E", "E": "S", "S": "W", "W": "N"},  # rot90
    {"N": "S", "E": "W", "S": "N", "W": "E"},  # rot180
    {"N": "W", "E": "N", "S": "E", "W": "S"},  # rot270
    # reflections (mirror x -> -x), then rotations
    {"N": "N", "E": "W", "S": "S", "W": "E"},  # mirror
    {"N": "W", "E": "S", "S": "E", "W": "N"},  # mirror+rot90
    {"N": "S", "E": "E", "S": "N", "W": "W"},  # mirror+rot180
    {"N": "E", "E": "N", "S": "W", "W": "S"},  # mirror+rot270
]


def all_cyclic_shifts(word: Sequence[str]) -> Iterable[List[str]]:
    w = len(word)
    for s in range(w):
        yield list(word[s:]) + list(word[:s])


def canonical_word(word: Sequence[str]) -> str:
    """
    Canonical representative under:
      - dihedral actions on directions (D4)
      - reversal+inversion
      - cyclic shifts
    Returned as a plain string of letters (expanded).
    """
    candidates: List[str] = []
    for mp in D4_MAPS:
        w1 = apply_dir_map(word, mp)
        w2 = reversed_inverted(w1)
        for wv in (w1, w2):
            for sh in all_cyclic_shifts(wv):
                candidates.append("".join(sh))
    return min(candidates)


# =============================================================================
# 7) Word scan (search) for a fixed torus
# =============================================================================

def scan_words(
    Lx: int,
    Ly: int,
    min_len: int = 4,
    max_len: int = 8,
    fix_first: str = "N",
    no_backtrack_flag: bool = True,
    symmetry_reduce: bool = True,
    require_distinct_Q: bool = True,
    require_full_support: bool = True,
    layout_mode: str = "row",
    dist_max: int = 3,
    min_k: int = 1,
    max_coset_layouts: int = 256,
    top_keep: int = 50,
) -> List[Dict[str, object]]:
    """
    Exhaustively enumerate words of lengths [min_len, max_len] and keep the best candidates.

    Returns a list of dict rows (ready to print or save).
    """
    if Lx % 2 or Ly % 2:
        raise ValueError("Require Lx, Ly even.")

    alphabet = list("NESW")
    seen = set()
    collected: List[CodeResult] = []

    for wlen in range(min_len, max_len + 1):
        tail_len = wlen - 1
        for tail in itertools.product(alphabet, repeat=tail_len):
            word = [fix_first] + list(tail)

            if no_backtrack_flag and not no_immediate_backtrack(word):
                continue

            if symmetry_reduce:
                can = canonical_word(word)
                if can in seen:
                    continue
            else:
                can = ""

            offsets = route_offsets(word)
            if require_distinct_Q and not offsets_distinct(offsets):
                continue
            if require_full_support and support_weight_mod2(offsets) != len(offsets):
                continue

            # Evaluate word (maybe multiple layouts)
            res_list = evaluate_word("".join(word), Lx, Ly, layout_mode=layout_mode, dist_max=dist_max, max_coset_layouts=max_coset_layouts)
            for res in res_list:
                if res.k >= min_k:
                    collected.append(res)

            if symmetry_reduce:
                # mark canonical word as seen once we've checked it
                seen.add(can)

    # convert results to rows + sort
    def score(res: CodeResult) -> Tuple[int, int, int]:
        # sort primarily by "lower bound" on distance
        if dist_max > 0:
            dx = dist_max + 1 if res.dX is None else res.dX
            dz = dist_max + 1 if res.dZ is None else res.dZ
            d = min(dx, dz)
        else:
            d = 0
        return (d, res.k, -res.word_len)

    collected.sort(key=score, reverse=True)

    rows: List[Dict[str, object]] = []
    for res in collected[:top_keep]:
        rows.append({
            "word": res.word,
            "w": res.word_len,
            "|P|": res.support_wt,
            "layout": res.layout,
            "Lx": res.Lx,
            "Ly": res.Ly,
            "n": res.n,
            "k": res.k,
            "dX": fmt_d(res.dX, dist_max),
            "dZ": fmt_d(res.dZ, dist_max),
        })

    return rows


# =============================================================================
# 8) Spyder entry points (edit CONFIG below and press F5)
# =============================================================================

def run_bench(
    words: List[str],
    tori: List[Tuple[int, int]],
    layout_mode: str = "row",
    dist_max: int = 0,
    max_rows_print: int = 60,
    save_csv: bool = False,
    csv_path: str = "out/bench_results.csv",
) -> None:
    """
    Evaluate chosen words on chosen tori; print a table.
    """
    all_rows: List[Dict[str, object]] = []
    for w in words:
        for (Lx, Ly) in tori:
            res_list = evaluate_word(w, Lx, Ly, layout_mode=layout_mode, dist_max=dist_max)
            for res in res_list:
                all_rows.append({
                    "word": res.word,
                    "w": res.word_len,
                    "|P|": res.support_wt,
                    "Lx": res.Lx,
                    "Ly": res.Ly,
                    "layout": res.layout,
                    "n": res.n,
                    "k": res.k,
                    "rankX": res.rankX,
                    "rankZ": res.rankZ,
                    "dX": fmt_d(res.dX, dist_max),
                    "dZ": fmt_d(res.dZ, dist_max),
                })

    # nice sort: by (k desc, d desc, w asc)
    def sort_key(r: Dict[str, object]):
        k = int(r["k"])
        # parse d
        def parse_d(x: object) -> int:
            s = str(x)
            if s.startswith(">"):
                return dist_max + 1
            if s == "-" or s == "":
                return 0
            return int(s)
        d = min(parse_d(r["dX"]), parse_d(r["dZ"]))
        return (-k, -d, int(r["w"]), str(r["word"]), int(r["Lx"]), int(r["Ly"]))

    all_rows.sort(key=sort_key)

    print_table(
        all_rows,
        columns=["word", "w", "|P|", "Lx", "Ly", "layout", "n", "k", "dX", "dZ", "rankX", "rankZ"],
        title=f"BENCH results (layout={layout_mode}, dist_max={dist_max})",
        max_rows=max_rows_print,
    )

    if save_csv:
        maybe_write_csv(csv_path, all_rows)


def run_search(
    Lx: int,
    Ly: int,
    min_len: int = 4,
    max_len: int = 8,
    fix_first: str = "N",
    no_backtrack_flag: bool = True,
    symmetry_reduce: bool = True,
    require_distinct_Q: bool = True,
    require_full_support: bool = True,
    layout_mode: str = "row",
    dist_max: int = 3,
    min_k: int = 1,
    top_keep: int = 50,
    max_rows_print: int = 30,
    save_csv: bool = False,
    csv_path: str = "out/search_results.csv",
) -> None:
    """
    Exhaustive scan over words and print top candidates.
    """
    rows = scan_words(
        Lx=Lx,
        Ly=Ly,
        min_len=min_len,
        max_len=max_len,
        fix_first=fix_first,
        no_backtrack_flag=no_backtrack_flag,
        symmetry_reduce=symmetry_reduce,
        require_distinct_Q=require_distinct_Q,
        require_full_support=require_full_support,
        layout_mode=layout_mode,
        dist_max=dist_max,
        min_k=min_k,
        top_keep=top_keep,
    )

    print_table(
        rows,
        columns=["word", "w", "|P|", "layout", "Lx", "Ly", "n", "k", "dX", "dZ"],
        title=f"SEARCH results on {Lx}x{Ly} (layout={layout_mode}, w={min_len}..{max_len}, dist_max={dist_max})",
        max_rows=max_rows_print,
    )

    if save_csv:
        maybe_write_csv(csv_path, rows)


# =============================================================================
# 9) EDIT THIS CONFIG IN SPYDER AND RUN (F5)
# =============================================================================

if __name__ == "__main__":

    # -------------------------
    # Choose a mode:
    #   MODE = "BENCH"  -> evaluate a list of words on a list of tori
    #   MODE = "SEARCH" -> exhaustive scan for new candidate words on one torus
    # -------------------------
    MODE = "BENCH"

    # -------------------------
    # Common parameters
    # -------------------------
    LAYOUT_MODE = "row"        # "row", "coset", or "row+coset"
    DIST_MAX = 4              # 0 disables distance search; otherwise exact up to this weight
    SAVE_CSV = False          # optional
    OUTDIR = "out"            # only used if SAVE_CSV=True

    # -------------------------
    # (A) BENCH configuration
    # -------------------------
    BENCH_WORDS = [
        "NE^2NE^2N",  # manuscript case study
        "NE3N",       # Riverlane benchmark
        "N2E3N2",     # Riverlane benchmark
    ]
    BENCH_TORI = [
        (12, 6),
        (16, 8),
        (24, 12),
    ]

    # -------------------------
    # (B) SEARCH configuration
    # -------------------------
    SEARCH_LX, SEARCH_LY = 16, 8
    SEARCH_MIN_LEN = 4
    SEARCH_MAX_LEN = 8
    SEARCH_TOP_KEEP = 50
    SEARCH_MIN_K = 1
    SEARCH_NO_BACKTRACK = True
    SEARCH_SYMMETRY_REDUCE = True
    SEARCH_DISTINCT_OFFSETS = True
    SEARCH_FULL_SUPPORT = True

    # -------------------------
    # Run
    # -------------------------
    if MODE.upper() == "BENCH":
        run_bench(
            words=BENCH_WORDS,
            tori=BENCH_TORI,
            layout_mode=LAYOUT_MODE,
            dist_max=DIST_MAX,
            max_rows_print=80,
            save_csv=SAVE_CSV,
            csv_path=os.path.join(OUTDIR, "bench_results.csv"),
        )
    elif MODE.upper() == "SEARCH":
        run_search(
            Lx=SEARCH_LX,
            Ly=SEARCH_LY,
            min_len=SEARCH_MIN_LEN,
            max_len=SEARCH_MAX_LEN,
            fix_first="N",
            no_backtrack_flag=SEARCH_NO_BACKTRACK,
            symmetry_reduce=SEARCH_SYMMETRY_REDUCE,
            require_distinct_Q=SEARCH_DISTINCT_OFFSETS,
            require_full_support=SEARCH_FULL_SUPPORT,
            layout_mode=LAYOUT_MODE,
            dist_max=DIST_MAX,
            min_k=SEARCH_MIN_K,
            top_keep=SEARCH_TOP_KEEP,
            max_rows_print=40,
            save_csv=SAVE_CSV,
            csv_path=os.path.join(OUTDIR, "search_results.csv"),
        )
    else:
        raise ValueError("MODE must be 'BENCH' or 'SEARCH'.")
