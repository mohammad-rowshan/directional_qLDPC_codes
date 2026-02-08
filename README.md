# Directional-code Numerics

This folder contains a Python script for **static** (algebraic) numerics for **directional (route-generated) CSS qLDPC codes** in the **square-grid checkerboard model**.


---

## 1) Requirements

- Python **3.9+** recommended (works best on 3.10/3.11)
- No external packages required (standard library only)

If you want, you *can* use pandas in Spyder, but the script does **not** depend on it.

---

## 2) Files

- `dircode_numerics.py`  

---

## 3) How to run 

1. Open the script file.
2. Scroll to the bottom section:

   ```python
   if __name__ == "__main__":
       MODE = "BENCH"  # or "SEARCH"
       ...
   ```

3. Edit the parameters you want.
4. Run file.


---

## 4) Two main modes

### A) BENCH mode (evaluate chosen words on chosen tori)

Set:

```python
MODE = "BENCH"
BENCH_WORDS = ["NE^2NE^2N", "NE3N", "N2E3N2"]
BENCH_TORI = [(12,6), (16,8), (24,12)]
LAYOUT_MODE = "row"      # or "coset" / "row+coset"
DIST_MAX = 4             # 0 disables distance search
```

What you get:
- A table with columns like `word`, `w`, `|P|`, `Lx`, `Ly`, `n`, `k`, `dX`, `dZ`, ...

Interpretation of `dX` / `dZ`:
- If `DIST_MAX = 0`, distances are not computed and shown as `-`
- If `DIST_MAX > 0`:
  - a number like `3` means the exact distance is `3` (found ≤ cutoff)
  - an entry like `>4` means **no nontrivial logical was found up to weight 4**, so the distance is at least 5

### B) SEARCH mode (scan for new candidate words on a fixed torus)

Set:

```python
MODE = "SEARCH"
SEARCH_LX, SEARCH_LY = 16, 8
SEARCH_MIN_LEN = 4
SEARCH_MAX_LEN = 8
DIST_MAX = 4
LAYOUT_MODE = "row"
```

The scan:
- enumerates words from length `SEARCH_MIN_LEN` to `SEARCH_MAX_LEN`
- fixes the first step to `N` (reduces redundancy)
- can reject immediate backtracking (`...NS...`, `...EW...`)
- can quotient by symmetry (dihedral + reversal/inversion + cyclic shifts)
- can require all offsets `Q_j` to be distinct and non-cancelling

You will see the top candidates by a simple score (higher distance lower bound first, then larger `k`, then shorter word).

---

## 5) Layout options

`LAYOUT_MODE` can be:

- `"row"`  
  The simple row-alternating layout used in the manuscript (X on even-y ancillas, Z on odd-y ancillas).

- `"coset"`  
  Enumerate layouts that are constant on cosets implied by the word’s **Δ_odd(W)** generators (layout-coset theorem idea).
  This can improve results but costs more compute time.

- `"row+coset"`  
  Evaluate both the row layout and all coset layouts.

**Important:** The script caps coset enumeration by `max_coset_layouts` (default 256).  
If the number of coset assignments is larger than that, it prints a warning and skips coset enumeration.

---

## 6) Practical tips (avoid slow runs)

- Start with `DIST_MAX = 3` or `4` for `n = 64`
- For `n > 150`, keep `DIST_MAX <= 3` (distance search is combinatorial)
- Use `SEARCH_MAX_LEN <= 8` initially
- Keep `LAYOUT_MODE="row"` for first scans; switch to `"row+coset"` later for the best candidates

---

## 7) Optional CSV export

At the bottom:

```python
SAVE_CSV = True
OUTDIR = "out"
```

This will save:
- `out/bench_results.csv` for BENCH mode
- `out/search_results.csv` for SEARCH mode
