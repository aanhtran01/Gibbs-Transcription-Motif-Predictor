# -*- coding: utf-8 -*-
import numpy as np
import random
import argparse

# Create the argument parser
parser = argparse.ArgumentParser(description='File paths')

# Add the file path arguments
parser.add_argument('bound_file', type=str, help='Path to the bound file')
parser.add_argument('test_file', type=str, help='Path to the test file')

# Parse the command-line arguments
args = parser.parse_args()

#function to read in fasta files 
def read_fasta_file(filename):
    with open(filename, 'r') as f:
        read_id = None
        read_seq = ""
        for line in f:
            if line.startswith(">"):
                if read_id is not None:
                    yield (read_id, read_seq)
                read_id = line.strip()[1:]
                read_seq = ""
            else:
                read_seq += line.strip()
        if read_id is not None:
            yield (read_id, read_seq)

#Computes the profile matrix of a given set of motifs.
def compute_profile_matrix(motifs):
    motif_length = len(motifs[0])  # Length of each motif
    num_motifs = len(motifs)  # Number of motifs
    motif_scale = 1 / num_motifs  # Frequency scale factor

    seq_to_numeric = 'ACGT0123'  # Mapping of nucleotides to numeric values
    numeric_dict = {seq_to_numeric[i]: int(seq_to_numeric[i + 4]) for i in range(4)}

    # Initialize the profile matrix with zeros
    profile_matrix = [[0 for _ in range(motif_length)] for __ in range(4)]

    # Count the occurrences of each nucleotide at each position in the motifs
    for motif in motifs:
        for i in range(motif_length):
            nucleotide = motif[i]
            numeric_value = numeric_dict[nucleotide]
            profile_matrix[numeric_value][i] += motif_scale

    return profile_matrix


#Computes the pseudocount profile matrix of a given set of motifs.
def compute_pseudocount_profile_matrix(motifs):
    motif_length = len(motifs[0])  # Length of each motif
    num_motifs = len(motifs)  # Number of motifs
    motif_scale = 1 / (num_motifs + 4)  # Frequency scale factor, accounting for pseudocounts

    seq_to_numeric = 'ACGT0123'  # Mapping of nucleotides to numeric values
    numeric_dict = {seq_to_numeric[i]: int(seq_to_numeric[i + 4]) for i in range(4)}

    # Initialize the pseudocount profile matrix with scaled pseudocounts
    pseudocount_matrix = [[motif_scale for _ in range(motif_length)] for __ in range(4)]

    # Count the occurrences of each nucleotide at each position in the motifs and apply pseudocounts
    for motif in motifs:
        for i in range(motif_length):
            nucleotide = motif[i]
            numeric_value = numeric_dict[nucleotide]
            pseudocount_matrix[numeric_value][i] += motif_scale

    return pseudocount_matrix


#Computes the score of a given set of motifs.
def compute_motif_score(motifs):
    motif_length = len(motifs[0])  # Length of each motif
    num_motifs = len(motifs)  # Number of motifs

    seq_to_numeric = 'ACGT0123'  # Mapping of nucleotides to numeric values
    numeric_dict = {seq_to_numeric[i]: int(seq_to_numeric[i + 4]) for i in range(4)}

    # Initialize the count matrix with zeros
    count_matrix = [[0 for _ in range(4)] for __ in range(motif_length)]

    # Count the occurrences of each nucleotide at each position in the motifs
    for motif in motifs:
        for i in range(motif_length):
            nucleotide = motif[i]
            numeric_value = numeric_dict[nucleotide]
            count_matrix[i][numeric_value] += 1

    score = 0
    # Compute the score as the sum of differences between the number of motifs and the maximum count at each position
    for i in range(motif_length):
        score += num_motifs - max(count_matrix[i])

    return num_motifs * motif_length - score


            
#Generates a random k-mer from a given sequence based on a profile matrix.
def generate_random_kmer_from_profile(sequence, motif_length, profile):
    sequence_length = len(sequence)  # Length of the sequence
    pattern_probabilities = []  # List to store pattern probabilities

    seq_to_numeric = 'ACGT0123'  # Mapping of nucleotides to numeric values
    seq_dict = {seq_to_numeric[i]: int(seq_to_numeric[i + 4]) for i in range(4)}

    # Compute the probabilities of patterns in the sequence based on the profile matrix
    for i in range(sequence_length - motif_length + 1):
        pattern = sequence[i:i + motif_length]
        pattern_probability = 1
        for j in range(motif_length):
            nucleotide = pattern[j]
            numeric_value = seq_dict[nucleotide]
            pattern_probability *= profile[numeric_value][j]
        pattern_probabilities.append(pattern_probability)

    total_probability = sum(pattern_probabilities)

    # Generate a random index from a biased probability distribution
    number = random.uniform(0, total_probability)
    cumulative_probability = 0
    for i, bias in enumerate(pattern_probabilities):
        cumulative_probability += bias
        if number <= cumulative_probability:
            return sequence[i:i + motif_length]


#Run Gibbs Sampling and return the best motifs and best score 
def run_gibbs_sampler(sequence_set, motif_length, num_sequences, num_iterations):
    sequence_length = len(sequence_set[0])  # Length of the sequences in the set
    motif_start_positions = [random.randint(0, sequence_length - motif_length) for _ in range(num_sequences)]  # Randomly select starting positions for motifs
    motifs = [sequence_set[i][motif_start_positions[i]:motif_start_positions[i] + motif_length] for i in range(num_sequences)]  # Extract motifs based on the starting positions
    best_motifs = motifs  # Initialize the best motifs with the initial motifs
    best_score = compute_motif_score(best_motifs)  # Compute the score of the best motifs using the motif score function

    for j in range(num_iterations):
        i = random.randint(0, num_sequences-1)  # Choose a random index i
        P = compute_pseudocount_profile_matrix(sequence_set[:i] + sequence_set[i+1:])  # Compute the pseudocount profile matrix without the i-th sequence
        motif_i = generate_random_kmer_from_profile(sequence_set[i], motif_length, P)  # Generate a random k-mer from the i-th sequence based on the profile matrix
        motifs = motifs[:i] + [motif_i] + motifs[i+1:]  # Update the motifs by replacing the i-th motif with the newly generated motif
        current_score = compute_motif_score(motifs)  # Compute the score of the updated motifs

        if current_score < best_score:
            best_motifs = motifs  # Update the best motifs if the current score is better
            best_score = current_score  # Update the best score

    return best_motifs, best_score



#Runs the Gibbs sampler multiple times and returns the best PWM for forward compliments.
def run_multiple_gibbs_sampler(dna_sequences, motif_length, num_sequences, num_iterations, _iter= 50):
    best_score = float('inf')  # Initialize the best score with infinity
    seed_value = 1234
    random.seed(seed_value)  # Seed the random number generator
    final_profile_matrix = None  # Initialize the final profile matrix

    for _ in range(_iter):
        curr_best_motifs, curr_best_score = run_gibbs_sampler(dna_sequences, motif_length, num_sequences, num_iterations)

        if curr_best_score < best_score:
            best_motifs, best_score = (curr_best_motifs, curr_best_score)  # Update the best motifs and score
            final_profile_matrix = compute_pseudocount_profile_matrix(best_motifs)  # Update the final profile matrix

    return final_profile_matrix


#function to find the reverse complement of a sequence
def ReverseComplement(seq):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement_seq = [complement_dict[base] for base in reversed(seq)]
    complement = ''.join(complement_seq)
    return complement

#function to calculate proabilities for kmers 
def calculate_probability(kmer, pwm):
    seq1 = 'ACGT0123'
    seq_dict = {seq1[i]: i for i in range(4)}
    p = 1
    for i in range(len(kmer)):
        nucleotide = kmer[i]
        index = seq_dict[nucleotide]
        p *= pwm[index][i]
    return p

#function to keep the highest probability among a kmer in a sequence 
def find_max_probability(sequence, kmer_length, pwm):
    max_prob = 0.0

    for i in range(len(sequence) - kmer_length + 1):
        kmer = sequence[i:i+kmer_length]
        probability = calculate_probability(kmer, pwm)
        if probability > max_prob:
            max_prob = probability

    return max_prob

#function to compute the highest probability for each forward sequence 
def all_forward_probabilites(test_reads, kmer_length, pwm):
    all_forward_probs = []
    
    for read_id, read_seq in test_reads:
        max_prob = find_max_probability(read_seq, kmer_length, pwm)
        read_probs = {'read_id': read_id, 'max_prob': max_prob}
        all_forward_probs.append(read_probs)
    
    return all_forward_probs


#function to compute the highest probability for each reverse sequence 
def all_reverse_probabilites(test_reads, kmer_length, pwm):
    all_rev_probs = []
    
    for read_id, read_seq in test_reads:
        rev_read = ReverseComplement(read_seq)
        max_prob = find_max_probability(rev_read, kmer_length, pwm)
        read_probs = {'read_id': read_id, 'max_prob': max_prob}
        all_rev_probs.append(read_probs)
    
    return all_rev_probs


# Function to merge and write the top 2000 sequences with the highest probabilities to a file
def merge_and_write_top_sequences(forward_probs, reverse_probs, output_file, top_k=2000):
    # Merge forward and reverse probabilities
    merged_probs = forward_probs + reverse_probs

    # Sort the sequences by max_prob in descending order
    merged_probs.sort(key=lambda x: x['max_prob'], reverse=True)

    # Write the top sequences to the output file
    written_seqs = set()  # Track sequences that have been written
    num_written_seqs = 0  # Track the number of written sequences

    with open(output_file, 'w') as file:
        i = 0
        while num_written_seqs < top_k and i < len(merged_probs):
            sequence = merged_probs[i]['read_id']
            if sequence not in written_seqs:
                file.write(f"{sequence}\n")
                written_seqs.add(sequence)
                num_written_seqs += 1
            i += 1


# Access the file paths using argparse
bound_file = args.bound_file
test_file = args.test_file


bound_reads = list(read_fasta_file(bound_file))
num_sequences = len(bound_reads)
#num_sequences = 500
#print(num_sequences)
motif_length = 16
num_iterations = 200

# Extract sequences from bound_reads
sequences = [seq for _, seq in bound_reads]

# Create a list of all reverse complements
all_sequences = sequences + [ReverseComplement(seq) for seq in sequences]

# Call run_multiple_gibbs_sampler with the sequences
PWM = run_multiple_gibbs_sampler(all_sequences, motif_length, len(all_sequences), num_iterations)

output_file = 'predictions.txt'

test_reads = read_fasta_file(test_file)
for_probs = all_forward_probabilites(test_reads, motif_length, PWM) 

test_reads = read_fasta_file(test_file)
rev_probs = all_reverse_probabilites(test_reads, motif_length, PWM) 

# Output a predictions.txt file of the top 2000 sequences that have the highest probability of binding to a specific transcription factor
merge_and_write_top_sequences(for_probs, rev_probs, output_file, top_k=2000)

