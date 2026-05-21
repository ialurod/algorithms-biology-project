import random

header = ">HBA_divergent_"
# Original Human HBA reference sequence
sequence = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
amino_acids = "ACDEFGHIKLMNPQRSTVWY"

# Target mutation rates: 10%, 20%, and 30%
rates = {"10": 0.10, "20": 0.20, "30": 0.30}

for level, rate in rates.items():
    seq_list = list(sequence)
    num_mutations = int(len(sequence) * rate)
    indices = random.sample(range(len(sequence)), num_mutations)
    
    for idx in indices:
        original = seq_list[idx]
        seq_list[idx] = random.choice([aa for aa in amino_acids if aa != original])
    
    # Export to FASTA format
    with open(f"HBA_divergent_{level}.fasta", "w") as f:
        f.write(f"{header}{level}\n")
        f.write("".join(seq_list) + "\n")