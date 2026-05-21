> **Use a simplified BLAST-like seed-and-extend implementation.**

We implemented a basic BLAST-like seed-and-extend algorithm based on the implementation covered in the practical lessons, with the following nuances:
* **Scoring:** Instead of a simple match/mismatch scoring system, we used the **BLOSUM62** substitution matrix to evaluate initial K-mer similarity. To avoid computational bottlenecks, we implemented a dictionary-based caching system for rapid score lookups.
* **Extension:** We used a strictly ungapped extension approach, keeping the drop-off threshold fixed at 5.
* We used Python’s built-in `time` module to compute runtimes.

---

> **Build a sequence database composed of proteins from a single family (e.g., retrieved from BLAST or Pfam) together with a small number of outliers.**

Our database was built by retrieving the first 100 entries that show up in UniProt when querying *"hemoglobin subunit alpha"* in the "Protein Name" advanced search. We downloaded these as a FASTA file and manually modified it to include five *hemoglobin subunit beta* and five *hemoglobin subunit delta* sequences from different species to act as our out-group outliers.

*The final database is accessible in this GitHub repository, along with R code and the algorithm.*

---

> **Systematically vary the k-mer size and seed threshold, and evaluate their impact on runtime and on the quality of the results obtained.**

We used human hemoglobin subunit alpha (sp|P69905) as the initial baseline query, alongside synthetically mutated versions at 10%, 20%, and 30% divergence. We stored the output in a CSV file and processed the results using several R packages (`ggplot2`, `dplyr`, etc.) to generate our visualizations.

<div align="center">
  <img width="1240" alt="Tradeoff Figure" src="plot_global_balance.png">
</div>

## Trade-off Analysis

After generating the similarity score data matrices and R visualizations, we observed the following parameter trade-offs:

* **Impact on Quality (Sensitivity vs. Specificity):** As the query diverges further from the database, the total number of valid `globin_hits` drops. A permissive BLOSUM seed threshold (T = 10) allows recovering highly mutated sequences (e.g., over 90 hits even at 30% divergence) but sacrifices specificity by capturing false positives (Beta/Delta outliers). Conversely, a strict threshold (T = 19) guarantees zero noise but prematurely rejects valid distant homologs.
* **Impact on Runtime (Efficiency):** Execution time is highly dependent on the number of seeds evaluated. Smaller K-mer sizes (K = 3) force the algorithm to evaluate vastly more initial sequence pairs, increasing the computational burden. However, our dictionary-based caching system for BLOSUM62 drastically reduced the exponential lookup overhead. Furthermore, stricter seed thresholds actively decrease runtime by filtering out weak matches before the extension phase is ever triggered.
* **Optimal Parameter Configuration:** To properly balance computational efficiency and alignment accuracy on divergent sequences, we conclude that the optimal parameter profile is **K-mer = 4** and **BLOSUM62 Seed Threshold = 16**. This configuration successfully recovers a high proportion of true globin hits across all divergence levels (safeguarding the detection of distant homologs) while strictly filtering out structural outliers (0 false positives) and maintaining highly stable and optimized execution times (~1.5s).
