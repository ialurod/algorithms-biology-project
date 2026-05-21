> **Use a simplified BLAST-like seed-and-extend implementation.**

We implemented a basic BLAST-like seed-and-extend algorithm based on the implementation covered in the practical lessons, with the following nuances:
* **Scoring:** Instead of using a BLOSUM62 matrix or a similar substitution matrix, we used a simple match/mismatch scoring system.
* **Extension:** We used a strictly ungapped extension approach.
* We used Python’s built-in `time` module to compute runtimes.

---

> **Build a sequence database composed of proteins from a single family (e.g., retrieved from BLAST or Pfam) together with a small number of outliers.**

Our database was built by retrieving the first 100 entries that show up in UniProt when querying *"hemoglobin subunit alpha"* in the "Protein Name" advanced search. We downloaded these as a FASTA file and manually modified it to include five *hemoglobin subunit beta* and five *hemoglobin subunit delta* sequences from different species to act as our out-group outliers.

*The final database is accessible in this GitHub repository, along with R code and the algorithm.*

---

> **Systematically vary the k-mer size and seed threshold, and evaluate their impact on runtime and on the quality of the results obtained.**

We used human hemoglobin subunit alpha as the initial query, stored the output in a CSV file, and processed the results using several R packages.

<div align="center">
  <img width="1240" alt="Tradeoff Figure" src="https://github.com/user-attachments/assets/cfedf03e-c0ff-4a23-a145-fa86cd0dd3b8" />
</div>

As captured in the figure above, there is a significant tradeoff between sensitivity and specificity. To capture all the target hemoglobin subunit alpha sequences, we must broaden our drop-off threshold, which inevitably results in capturing an outlier. However, this is an acceptable compromise because we recover a significant number of true positive hits in exchange. 

Factoring in the computational runtimes, we concluded that for this specific case, **a drop-off threshold of 3 and a k-mer size of 3 is optimal.**

### Investigating False Negatives
During our analysis, we noticed that 5 target HBA sequences were consistently escaping capture. We investigated the UniProt IDs and identified the following species:

* `sp|P02020|HBA_LEPPA`
* `sp|P02021|HBA_HETPO`
* `sp|P56691|HBA_HEMAK`
* `sp|P84216|HBA_BATEA`
* `sp|Q9YGW2|HBA_MUSGR`

These sequences correspond to two distinct evolutionary extremes:
1. **Three shark species and a lungfish:** These ancient lineages diverged from the human common ancestor roughly 450 million years ago, resulting in massive amounts of mutations and sequence indels (insertions/deletions). Because our algorithm relies on *ungapped* extension, it cannot successfully traverse these indels. Consequently, the drop-off threshold is triggered, and the algorithm terminates prematurely.
2. **An extremophile dragonfish:** The HBA of this species has been subjected to extreme natural selection over millions of years to function in sub-zero temperatures, leading to wild structural sequence differences when compared to human HBA.

Taking into account the structural limitations of our simplified BLAST implementation, we can confidently conclude that the algorithm is correctly demonstrating expected biological limitations and is working as intended.

### Investigating False Positive

HBB_TREBE, another fish hemoglobin beta subunit, is consistently passed with dropoff threshold > 2 because it contains a highly conserved paralog domain with human HBA that strictly meets the dropoff threshold everytime. It is then, a loose outlier.

> **Perform the analysis on query sequences with different levels of divergence and determine parameter ranges that balance computational efficiency and alignment accuracy.**

![Divergence Analysis Performance](plot_divergence.png)

To systematically evaluate the algorithm's performance under evolutionary stress, we generated three divergent queries based on the human HBA sequence, introducing random point mutations at rates of 10%, 20%, and 30%. Based on the generated data matrices and R visualizations, we observed the following parameter trade-offs:

* **Impact on Quality (Sensitivity):** As the query diverges further from the database, the total number of valid `globin_hits` drops significantly (from 95 in the baseline down to a maximum of 55 at 30% divergence). Strict drop-off thresholds (T = 2) prematurely terminate the ungapped extension phase when encountering localized mutation noise, causing further hit loss. A higher threshold (T >= 4) is necessary to tolerate these introduced mismatches and recover the maximum possible alignments.
* **Impact on Runtime (Efficiency):** Execution time decreases as divergence increases due to a severe reduction in exact initial K-mer matches. Furthermore, smaller K-mer sizes (K = 2) consistently bottleneck the algorithm due to the massive generation of initial seeds. For example, at 10% divergence with T = 2, increasing K from 2 to 5 cuts the runtime in half (from 0.0664s to 0.0302s). 
* **Optimal Parameter Configuration:** To properly balance computational efficiency and alignment accuracy on divergent sequences, we conclude that the optimal parameter range is **K-mer = 4** and **Drop-off Threshold = 4 or 5**. This configuration successfully recovers the maximum possible true globin hits across all divergence levels (e.g., stabilizing at 55 hits for the 30% divergent query) while avoiding the exponential lookup overhead of smaller seeds, maintaining highly optimized execution times (~0.02s).
