> **Use a simplified BLAST-like seed-and-extend implementation with a protein substitution matrix.**

We evolved our simplified BLAST-like seed-and-extend algorithm into a true similarity-driven heuristic by integrating the **BLOSUM62** substitution matrix directly into the seeding phase:
* **Similarity Seeding:** Instead of restricting the algorithm to exact K-mer matches, we compute the cumulative substitution score for every query and target K-mer pair using BLOSUM62. A seed is only considered a valid hit and passed to the extension phase if its total score satisfies a minimum **Seed Threshold (T)**.
* **Algorithmic Optimization:** To eliminate the severe performance bottlenecks caused by millions of matrix lookups within nested loops, we extracted the BLOSUM62 data into a native Python dictionary (`fast_blosum`) and built a **lazy-evaluation score cache**. This optimization completely bypassed redundant calculations, reducing execution time from minutes to seconds.
* **Extension:** We maintained a strictly ungapped extension approach, keeping the drop-off threshold fixed at 5 while shifting our experimental focus to the optimization of the seed threshold.

---

> **Build a sequence database composed of proteins from a single family (e.g., retrieved from BLAST or Pfam) together with a small number of outliers.**

Our database was built by retrieving the first 100 entries that show up in UniProt when querying *"hemoglobin subunit alpha"* in the "Protein Name" advanced search. We downloaded these as a FASTA file and manually modified it to include five *hemoglobin subunit beta* and five *hemoglobin subunit delta* sequences from different species to act as our out-group outliers.

*The final database is accessible in this GitHub repository, along with R code and the algorithm.*

---

> **Systematically vary the k-mer size and seed threshold, and evaluate their impact on runtime and on the quality of the results obtained.**

We used human hemoglobin subunit alpha (sp|P69905) as the initial baseline query (0% divergence), alongside synthetically mutated versions at 10%, 20%, and 30% sequence divergence. We systematically tested K-mer sizes ($K \in \{3, 4, 5\}$) against a range of BLOSUM62 Seed Thresholds ($T \in \{10, 13, 16, 19\}$). We stored the output in a CSV file and processed the results using R to analyze accuracy and computational efficiency.

<div align="center">
  <img width="1240" alt="Tradeoff Figure" src="plot_global_balance.png">
</div>

## Trade-off Analysis (BLOSUM62 Parameter Balancing)

After generating the empirical data matrices and R visualizations, we observed the following parameter trade-offs:

* **Sensitivity vs. Specificity (Quality):** A permissive seed threshold ($T = 10$) maximizes sensitivity by capturing distant homologs even at 30% query divergence (recovering up to 90 valid globin hits), but it suffers from poor specificity, allowing structural outliers (Beta/Delta hemoglobin "trampas") to break through as false positives. Conversely, a strict threshold ($T = 19$) removes all background noise but prematurely rejects true positive alignments on highly mutated queries.
* **Computational Efficiency (Runtime):** Execution times are strongly tied to the interaction between $K$ and $T$. Smaller K-mers ($K = 3$) generate an immense pool of initial seed candidates, forcing more frequent extensions. However, stricter seed thresholds ($T \ge 16$) act as an aggressive early-stage filter, killing weak alignments letter-by-letter before full extension is triggered, significantly optimizing total runtime.
* **Optimal Parameter Profile:** To achieve an ideal balance between sensitivity on divergent evolutionary lineages and algorithmic speed, we conclude that the optimal parameter profile is **K-mer = 4** and **BLOSUM62 Seed Threshold = 16**. This precise configuration filters out structural noise completely (0 false positives), maintains highly optimized search speeds (~1.5s), and preserves a robust capacity to recover true distant globin relationships across all divergence levels.
