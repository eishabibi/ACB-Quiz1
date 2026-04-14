import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# Q1: Approximate Pattern Matching (Hamming Distance d=1)
# ============================================================

def hamming_distance(s1, s2):
    """Calculate Hamming distance between two strings of equal length."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def approximate_pattern_matching(text, pattern, d):
    """
    Find all starting positions where pattern appears in text
    with at most d mismatches (Hamming distance <= d).
    """
    positions = []
    n = len(text)
    k = len(pattern)
    for i in range(n - k + 1):
        substring = text[i:i+k]
        if hamming_distance(substring, pattern) <= d:
            positions.append(i)
    return positions

print("=" * 60)
print("Q1: Approximate Pattern Matching")
print("=" * 60)
text    = "ATGCATATGACTACTAGATACTGATACTGATACATA"
pattern = "ATAG"
d       = 1

positions = approximate_pattern_matching(text, pattern, d)

print(f"Text    : {text}")
print(f"Pattern : {pattern}")
print(f"Hamming Distance (d): {d}")
print(f"\nMatching positions (0-indexed): {positions}")
print("\nVerification:")
for pos in positions:
    substr = text[pos:pos+len(pattern)]
    hd = hamming_distance(substr, pattern)
    print(f"  Position {pos}: '{substr}' -> Hamming distance = {hd}")


# ============================================================
# Q2: Cumulative GC Skew
# ============================================================

def generate_skewed_genome(length=50000, split=25000):
    np.random.seed(42)  # for reproducibility
    p1 = [0.45, 0.2, 0.3, 0.05]   # A, C, G, T  (C > G, Lagging strand)
    seq1 = np.random.choice(['A', 'C', 'G', 'T'], size=split, p=p1)
    p2 = [0.05, 0.3, 0.2, 0.45]   # A, C, G, T  (G > C, Leading strand)
    seq2 = np.random.choice(['A', 'C', 'G', 'T'], size=length-split, p=p2)
    return "".join(np.concatenate([seq1, seq2]))

def compute_gc_skew(genome):
    """Compute cumulative GC skew: skew[i] = #G - #C up to position i."""
    skew = [0]
    for nucleotide in genome:
        if nucleotide == 'G':
            skew.append(skew[-1] + 1)
        elif nucleotide == 'C':
            skew.append(skew[-1] - 1)
        else:
            skew.append(skew[-1])
    return skew

print("\n" + "=" * 60)
print("Q2: Cumulative GC Skew")
print("=" * 60)

genome = generate_skewed_genome()
print(f"Genome length     : {len(genome)} bp")
print(f"G count           : {genome.count('G')}")
print(f"C count           : {genome.count('C')}")

skew = compute_gc_skew(genome)
min_skew_pos = skew.index(min(skew))
max_skew_pos = skew.index(max(skew))

print(f"Min GC Skew value : {min(skew)} at position {min_skew_pos}")
print(f"Max GC Skew value : {max(skew)} at position {max_skew_pos}")
print(f"(Ori of replication likely near position {min_skew_pos})")

# Plot
plt.figure(figsize=(12, 5))
plt.plot(skew, color='steelblue', linewidth=0.8, label='Cumulative GC Skew')
plt.axvline(x=25000, color='red', linestyle='--', linewidth=1.2, label='Midpoint (bias shift)')
plt.axvline(x=min_skew_pos, color='orange', linestyle='--', linewidth=1.2,
            label=f'Min Skew (pos={min_skew_pos})')
plt.axvline(x=max_skew_pos, color='green', linestyle='--', linewidth=1.2,
            label=f'Max Skew (pos={max_skew_pos})')
plt.xlabel('Position in Genome (bp)')
plt.ylabel('Cumulative GC Skew (G - C)')
plt.title('Cumulative GC Skew of Synthetic Genome (50,000 bp)')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig('/home/claude/gc_skew_plot.png', dpi=150)
plt.close()
print("Plot saved as gc_skew_plot.png")


# ============================================================
# Q3: Motif Matrix Analysis
# ============================================================

print("\n" + "=" * 60)
print("Q3: Motif Matrix Analysis")
print("=" * 60)

motifs = [
    'AACCTA',
    'CCCGTT',
    'GACCTT',
    'GGATTA',
    'TTCCGG'
]

print(f"Motif Matrix:")
for m in motifs:
    print(f"  {m}")

k = len(motifs[0])
t = len(motifs)

# --- Count Matrix ---
def build_count_matrix(motifs):
    k = len(motifs[0])
    count = {n: [0]*k for n in 'ACGT'}
    for motif in motifs:
        for j, nuc in enumerate(motif):
            count[nuc][j] += 1
    return count

# --- Profile Matrix ---
def build_profile(motifs, pseudocount=0):
    t = len(motifs)
    count = build_count_matrix(motifs)
    profile = {}
    for nuc in 'ACGT':
        profile[nuc] = [(c + pseudocount) / (t + 4*pseudocount) for c in count[nuc]]
    return profile

# --- Consensus Sequence ---
def consensus_sequence(motifs):
    count = build_count_matrix(motifs)
    consensus = ''
    for j in range(len(motifs[0])):
        best_nuc = max('ACGT', key=lambda n: count[n][j])
        consensus += best_nuc
    return consensus

# --- Motif Score ---
def motif_score(motifs):
    count = build_count_matrix(motifs)
    t = len(motifs)
    score = 0
    for j in range(len(motifs[0])):
        max_count = max(count[n][j] for n in 'ACGT')
        score += (t - max_count)
    return score

# --- Probability of a k-mer ---
def kmer_probability(kmer, profile):
    prob = 1.0
    for j, nuc in enumerate(kmer):
        prob *= profile[nuc][j]
    return prob

# --- Most Probable k-mer ---
def most_probable_kmer(text, k, profile):
    best_kmer = text[:k]
    best_prob = -1
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        prob = kmer_probability(kmer, profile)
        if prob > best_prob:
            best_prob = prob
            best_kmer = kmer
    return best_kmer, best_prob

# --- Greedy Motif Search ---
def greedy_motif_search(dna, k, t):
    best_motifs = [seq[:k] for seq in dna]
    first_seq = dna[0]
    for i in range(len(first_seq) - k + 1):
        motifs = [first_seq[i:i+k]]
        for j in range(1, t):
            profile = build_profile(motifs, pseudocount=1)
            motifs.append(most_probable_kmer(dna[j], k, profile)[0])
        if motif_score(motifs) < motif_score(best_motifs):
            best_motifs = motifs
    return best_motifs

# ---- ANSWERS ----
consensus = consensus_sequence(motifs)
score     = motif_score(motifs)
profile   = build_profile(motifs, pseudocount=0)
count_mat = build_count_matrix(motifs)

print(f"\n--- Count Matrix ---")
print(f"{'':4}", end="")
for j in range(k):
    print(f"  Col{j+1}", end="")
print()
for nuc in 'ACGT':
    print(f"  {nuc} :", end="")
    for val in count_mat[nuc]:
        print(f"  {val:4}", end="")
    print()

print(f"\n--- Profile Matrix ---")
print(f"{'':4}", end="")
for j in range(k):
    print(f"  Col{j+1}", end="")
print()
for nuc in 'ACGT':
    print(f"  {nuc} :", end="")
    for val in profile[nuc]:
        print(f"  {val:.2f}", end="")
    print()

print(f"\n3.1 Consensus Sequence  : {consensus}")
print(f"3.2 Motif Score         : {score}")

profile_with_pseudo = build_profile(motifs, pseudocount=1)
print(f"\n3.3 Most Probable {k}-mer (per motif row using profile):")
for motif in motifs:
    kmer, prob = most_probable_kmer(motif, k, profile_with_pseudo)
    print(f"    In '{motif}': most probable k-mer = '{kmer}' (prob={prob:.5f})")

best_motifs   = greedy_motif_search(motifs, k, t)
best_score    = motif_score(best_motifs)
best_consensus = consensus_sequence(best_motifs)

print(f"\n3.4 Best Motifs (GreedyMotifSearch with Laplace pseudocounts):")
for bm in best_motifs:
    print(f"    {bm}")
print(f"    Score     : {best_score}")
print(f"    Consensus : {best_consensus}")

print("\n" + "=" * 60)
print("ALL DONE!")
print("=" * 60)
