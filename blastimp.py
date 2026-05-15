import time
from Bio import SeqIO

def load_database(filename):
    database = {}
    for record in SeqIO.parse(filename, "fasta"):
        database[record.id] = str(record.seq)
    return database

def get_kmers_with_positions(sequence, k):
    kmers = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer not in kmers:
            kmers[kmer] = []
        kmers[kmer].append(i)
    return kmers

def extend_seed(query_seq, target_seq, q_start, t_start, k, threshold):
    max_score = k
    current_score = k
    match_reward = 1
    mismatch_penalty = -1
    
    q_right = q_start + k
    t_right = t_start + k
    
    while q_right < len(query_seq) and t_right < len(target_seq):
        if query_seq[q_right] == target_seq[t_right]:
            current_score += match_reward
            max_score = max(max_score, current_score)
        else:
            current_score += mismatch_penalty
            
        if max_score - current_score >= threshold:
            break
            
        q_right += 1
        t_right += 1

    q_left = q_start - 1
    t_left = t_start - 1
    current_score = max_score
    
    while q_left >= 0 and t_left >= 0:
        if query_seq[q_left] == target_seq[t_left]:
            current_score += match_reward
            max_score = max(max_score, current_score)
        else:
            current_score += mismatch_penalty
            
        if max_score - current_score >= threshold:
            break
            
        q_left -= 1
        t_left -= 1

    return max_score

def run_blast_search(query_seq, database, k_size, threshold):
    start_time = time.time()
    query_seeds = get_kmers_with_positions(query_seq, k_size)
    results = {}
    
    for target_id, target_seq in database.items():
        best_target_score = 0
        target_kmers = get_kmers_with_positions(target_seq, k_size)
        
        for t_kmer, t_positions in target_kmers.items():
            if t_kmer in query_seeds:
                q_positions = query_seeds[t_kmer]
                for q_pos in q_positions:
                    for t_pos in t_positions:
                        score = extend_seed(query_seq, target_seq, q_pos, t_pos, k_size, threshold)
                        if score > best_target_score:
                            best_target_score = score
                            
        if best_target_score >= 10: 
            results[target_id] = best_target_score
            
    end_time = time.time()
    return results, end_time - start_time

if __name__ == "__main__":
    db = load_database("C:/Users/Ian/Documents/AB/Group Project/my_database.fasta")
    
    query_id = "sp|P69905|HBA_HUMAN"
    query_sequence = db[query_id]
    
    k_mer_sizes = [2, 3, 4, 5]
    thresholds = [2, 3, 4, 5] 

    with open("C:/Users/Ian/Documents/AB/Group Project/output_file.csv", "wt") as file:
        file.write(f"query_id, test_k, test_threshold, globin_hits, outlier_hits, runtime_seconds\n")
        for test_k in k_mer_sizes:
            for test_thresh in thresholds:
                hits, search_time = run_blast_search(query_sequence, db, test_k, test_thresh)
                globin_hits = sum(1 for seq_id in hits.keys() if "HBA" in seq_id)
                outlier_hits = sum(1 for seq_id in hits.keys() if "HBA" not in seq_id)
                for seq_id in hits.keys():
                    if "HBA" not in seq_id:
                        print("I am a captured outlier:", seq_id)

                all_alphas_in_db = [seq_id for seq_id in db.keys() if "HBA" in seq_id]

                missing_alphas = [seq_id for seq_id in all_alphas_in_db if seq_id not in hits]

                # 3. Reveal them!
                print("The 5 Rejected Alpha Globins are:")
                for missing in missing_alphas:
                    print(missing)
                file.write(f"{query_id}, {test_k}, {test_thresh}, {globin_hits}, {outlier_hits}, {search_time:.4f}\n")