Gibbs Transcription Motif Predictor
===================================

This is a Python script that predicts sequence motifs using Gibbs sampling and computes the highest probability of binding for a set of test sequences. The script takes input files containing bound sequences and test sequences in FASTA format.

Requirements
------------

Python 3.x
NumPy library
Installation

Make sure you have Python 3.x installed on your system.
Install the NumPy library by running the following command:
pip install numpy


Usage
-----

Import the required libraries at the beginning of your Python script:

import numpy as np
import random
import argparse

Create an instance of the argument parser and add the file path arguments:

parser = argparse.ArgumentParser(description='File paths')
parser.add_argument('bound_file', type=str, help='Path to the bound file')
parser.add_argument('test_file', type=str, help='Path to the test file')
args = parser.parse_args()

Function Definitions:
----------------------------------------------
read_fasta_file(filename): Reads sequences from a FASTA file and yields (read_id, read_seq) pairs.
compute_profile_matrix(motifs): Computes the profile matrix of a given set of motifs.
compute_pseudocount_profile_matrix(motifs): Computes the pseudocount profile matrix of a given set of motifs.
compute_motif_score(motifs): Computes the score of a given set of motifs.
generate_random_kmer_from_profile(sequence, motif_length, profile): Generates a random k-mer from a given sequence based on a profile matrix.
run_gibbs_sampler(sequence_set, motif_length, num_sequences, num_iterations): Runs the Gibbs sampler and returns the best motifs and best score.
run_multiple_gibbs_sampler(dna_sequences, motif_length, num_sequences, num_iterations, _iter=50): Runs the Gibbs sampler multiple times and returns the best PWM for forward compliments.
ReverseComplement(seq): Finds the reverse complement of a DNA sequence.
calculate_probability(kmer, pwm): Calculates the probability of a k-mer given a PWM.
find_max_probability(sequence, kmer_length, pwm): Finds the highest probability among k-mers in a sequence.
all_forward_probabilites(test_reads, kmer_length, pwm): Computes the highest probability for each forward sequence.
all_reverse_probabilites(test_reads, kmer_length, pwm): Computes the highest probability for each reverse sequence.
merge_and_write_top_sequences(forward_probs, reverse_probs, output_file, top_k=2000): Merges and writes the top sequences with the highest probabilities to a file.

Access the file paths using argparse and read the bound sequences:
------------------------------------------------------------------
bound_file = args.bound_file
test_file = args.test_file
bound_reads = list(read_fasta_file(bound_file))

Set the desired parameters such as motif length, number of iterations, and top sequences to write:
--------------------------------------------------------------------------------------------------
motif_length = 16
num_iterations = 200
output_file = 'predictions.txt'

Run the motif prediction and probability calculation:
-----------------------------------------------------
sequences = [seq for _, seq in bound_reads]
all_sequences = sequences + [ReverseComplement(seq) for seq in sequences]
PWM = run_multiple_gibbs_sampler(all_sequences, motif_length, len(all_sequences), num_iterations)
test_reads = read_fasta_file(test_file)
forward_probs = all_forward_probabilites(test_reads, motif_length, PWM)
test_reads = read_fasta_file(test_file)
reverse_probs = all_reverse_probabilites(test_reads, motif_length, PWM)

Merge the forward and reverse probabilities and write the top sequences to the output file:
-------------------------------------------------------------------------------------------
merge_and_write_top_sequences(forward_probs, reverse_probs, output_file, top_k=2000)
Example Usage

Assuming you have a bound file named "bound_sequences.fasta" and a test file named "test_sequences.fasta", you can run the script with the following command:
-------------------------------------------------------------------------------------------------------------------------------------------------------------
motif_predictor.py bound_sequences.fasta test_sequences.fasta

The script will generate a file named "predictions.txt" containing the top 2000 sequences with the highest probability of binding to a specific transcription factor.

Future Improvements:
--------------------
Due to the nature of Gibbs sampling using random starts, the parameters and interation size have not been tuned to produce optimal results. A larger number of iterations must be done in order to make sure that the PWM has converged. Program can be updated so that it is not dependent on iterations but rather on PWN convergence status.

Other methods such as an Expectation-Maximization or a Logistic Regression machine learning model might better accurately predict transcription motifs in comparison to the Gibbs Sampling methid. 