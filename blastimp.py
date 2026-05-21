import time
from Bio import SeqIO
from Bio.Align import substitution_matrices

matrix = substitution_matrices.load("BLOSUM62")
fast_blosum = {}
for (a, b), score in matrix.items():
    fast_blosum[(a, b)] = score
    fast_blosum[(b, a)] = score

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

def run_blast_search(query_seq, database, k_size, dropoff_thresh, seed_thresh):
    start_time = time.time()
    query_seeds = get_kmers_with_positions(query_seq, k_size)
    results = {}
    
    score_cache = {}
    
    for target_id, target_seq in database.items():
        best_target_score = 0
        target_kmers = get_kmers_with_positions(target_seq, k_size)
        
        for t_kmer, t_positions in target_kmers.items():
            for q_kmer, q_positions in query_seeds.items():
                
                pair_key = (q_kmer, t_kmer)
                if pair_key in score_cache:
                    kmer_score = score_cache[pair_key]
                else:
                    kmer_score = 0
                    for a, b in zip(q_kmer, t_kmer):
                        kmer_score += fast_blosum.get((a, b), -4)
                    score_cache[pair_key] = kmer_score
                
                if kmer_score >= seed_thresh:
                    for q_pos in q_positions:
                        for t_pos in t_positions:
                            score = extend_seed(query_seq, target_seq, q_pos, t_pos, k_size, dropoff_thresh)
                            if score > best_target_score:
                                best_target_score = score
                                
        if best_target_score >= 10: 
            results[target_id] = best_target_score
            
    end_time = time.time()
    return results, end_time - start_time

if __name__ == "__main__":
    db = load_database("my_database.fasta")
    
    with open("HBA_divergent_30.fasta", "r") as f:
        query_id = f.readline().strip()[1:] 
        query_sequence = f.readline().strip()
    
    k_mer_sizes = [3, 4, 5]
    dropoff_fixed = 5
    seed_thresholds = [10, 13, 16, 19]

    with open("output_file.csv", "a") as file:
        for test_k in k_mer_sizes:
            for test_seed in seed_thresholds:
                hits, search_time = run_blast_search(query_sequence, db, test_k, dropoff_fixed, test_seed)
                
                globin_hits = sum(1 for seq_id in hits.keys() if "HBA" in seq_id)
                outlier_hits = sum(1 for seq_id in hits.keys() if "HBA" not in seq_id)
                
                for seq_id in hits.keys():
                    if "HBA" not in seq_id:
                        print("I am a captured outlier:", seq_id)

                all_alphas_in_db = [seq_id for seq_id in db.keys() if "HBA" in seq_id]
                missing_alphas = [seq_id for seq_id in all_alphas_in_db if seq_id not in hits]

                print("The rejected Alpha Globins are:")
                for missing in missing_alphas:
                    print(missing)
                
                file.write(f"{query_id}, {test_k}, {test_seed}, {globin_hits}, {outlier_hits}, {search_time:.4f}\n")