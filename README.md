> **Use a simplified BLAST-like seed-and-extend implementation with a protein substitution matrix.**

We evolved our simplified BLAST-like seed-and-extend algorithm into a true similarity-driven heuristic by integrating the **BLOSUM62** substitution matrix directly into the seeding phase. To optimize performance, we implemented a lazy-evaluation cache using Python dictionaries, which reduced redundant matrix lookups and significantly improved execution speed.

---

> **Build a sequence database composed of proteins from a single family.**

Our database was built using UniProt entries for *"hemoglobin subunit alpha"*, augmented with five *hemoglobin subunit beta* and five *hemoglobin subunit delta* sequences to serve as structural outliers.

---

## Performance Analysis

To validate our implementation, we conducted a systematic analysis across three key dimensions:

### 1. Specificity & Baseline Noise Filtering
This analysis focuses on our control set (Human Alpha Hemoglobin) and the ability to exclude structural outliers (Beta/Delta chains). The chart below shows that low seed thresholds (T=10) allow significant background noise, whereas higher thresholds (T >= 16) successfully filter out false positives while maintaining high sensitivity.

![Baseline Performance](plot_baseline.png)

### 2. Sensitivity on Divergent Sequences
To test evolutionary robustness, we introduced synthetic mutations (10%, 20%, 30% divergence). This plot demonstrates the "evolutionary reach" of our algorithm. By adjusting the BLOSUM62 seed threshold, we can recover distant homologs that would be missed by exact-match algorithms.

![Divergence Analysis](plot_divergence.png)

### 3. Optimization & Global Trade-off
Finally, we analyzed the intersection of execution time (efficiency) and retrieval accuracy (sensitivity). This plot identifies the "sweet spot" for our parameters, showing how we balance computational costs with biological discovery.

![Global Balance](plot_global_balance.png)

**Conclusion:** Based on these results, we determined that **K-mer = 4** and **Seed Threshold = 16** provides the optimal balance, ensuring zero false positives, stable performance, and broad sequence coverage.
