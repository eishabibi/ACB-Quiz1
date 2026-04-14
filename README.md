# ACB Quiz 1 - Bioinformatics Solutions

Solutions to all three questions from Quiz 1 covering:
- Approximate Pattern Matching
- GC Skew Analysis
- Motif Matrix Analysis

---

## Q1: Approximate Pattern Matching

**Problem:** Find positions of Pattern `ATAG` in Text using Hamming distance d=1.

**Algorithm:** Slide a window of length `k` over the text. At each position, compute the Hamming distance (number of mismatches) between the window and the pattern. If distance ≤ d, record the position.

**Result:**
```
Text    : ATGCATATGACTACTAGATACTGATACTGATACATA
Pattern : ATAG
d       : 1

Matching positions (0-indexed): [4, 13, 17, 23, 29]

Verification:
  Position  4: 'ATAT' -> Hamming distance = 1
  Position 13: 'CTAG' -> Hamming distance = 1
  Position 17: 'ATAC' -> Hamming distance = 1
  Position 23: 'ATAC' -> Hamming distance = 1
  Position 29: 'ATAC' -> Hamming distance = 1
```

---

## Q2: Cumulative GC Skew

**Problem:** Generate a synthetic 50,000 bp genome with a bias shift at position 25,000 and plot the cumulative GC Skew.

**Algorithm:**
- GC Skew at position i: +1 if G, -1 if C, 0 otherwise
- Cumulative skew is the running sum
- Minimum of skew → likely origin of replication (ori)

**Result:**
```
Genome length     : 50,000 bp
Min GC Skew value : -324 at position 49957
Max GC Skew value :  2459 at position 24998
```

The plot (`gc_skew_plot.png`) clearly shows the skew rising in the first half (C > G region) and falling in the second half (G > C region), with the inflection near the midpoint.

---

## Q3: Motif Matrix Analysis

**Motif Matrix:**
```
AACCTA
CCCGTT
GACCTT
GGATTA
TTCCGG
```

### Count Matrix
| | Col1 | Col2 | Col3 | Col4 | Col5 | Col6 |
|---|---|---|---|---|---|---|
| A | 1 | 2 | 1 | 0 | 0 | 2 |
| C | 1 | 1 | 4 | 3 | 0 | 0 |
| G | 2 | 1 | 0 | 1 | 1 | 1 |
| T | 1 | 1 | 0 | 1 | 4 | 2 |

### Profile Matrix
| | Col1 | Col2 | Col3 | Col4 | Col5 | Col6 |
|---|---|---|---|---|---|---|
| A | 0.20 | 0.40 | 0.20 | 0.00 | 0.00 | 0.40 |
| C | 0.20 | 0.20 | 0.80 | 0.60 | 0.00 | 0.00 |
| G | 0.40 | 0.20 | 0.00 | 0.20 | 0.20 | 0.20 |
| T | 0.20 | 0.20 | 0.00 | 0.20 | 0.80 | 0.40 |

### 3.1 Consensus Sequence
```
GACCTA
```
(Most frequent nucleotide at each column)

### 3.2 Motif Score
```
Score = 13
```
(Total number of mismatches across all motifs vs consensus)

### 3.3 Most Probable k-mer
Using profile with Laplace pseudocounts (to avoid zero probabilities):
```
In 'AACCTA': AACCTA  (prob=0.00339)
In 'CCCGTT': CCCGTT  (prob=0.00113)
In 'GACCTT': GACCTT  (prob=0.00508)
In 'GGATTA': GGATTA  (prob=0.00068)
In 'TTCCGG': TTCCGG  (prob=0.00060)
```

### 3.4 Best Motifs (GreedyMotifSearch)
```
AACCTA
CCCGTT
GACCTT
GGATTA
TTCCGG
Score     : 13
Consensus : GACCTA
```

---

## Requirements

```
numpy
matplotlib
```

Install with:
```bash
pip install numpy matplotlib
```

## Run

```bash
python3 quiz_solutions.py
```
