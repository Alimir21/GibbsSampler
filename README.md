# GibbsSampler: A Motif Finding Algorithm

This project implements the GibbsSampler, a bioinformatics algorithm for finding motifs in DNA sequences. The algorithm iteratively refines a set of candidate motifs to find the most probable ones.

## Overview

The GibbsSampler is a probabilistic method that uses iterative sampling to improve motif selection. Motif finding is crucial in bioinformatics for identifying recurring patterns in biological sequences, such as protein binding sites.

## Functions

- `GibbsSampler(Dna, k, t, N)`: Main function to perform the Gibbs Sampler algorithm.
- `RandomMotifs(Dna, k, t)`: Generates random k-mers from each DNA string.
- `CountWithPseudocounts(Motifs)`: Counts occurrences of each nucleotide in the motifs, adding pseudocounts for stability.
- `ProfileWithPseudocounts(Motifs)`: Generates a profile matrix from motifs, incorporating pseudocounts.
- `Consensus(Motifs)`: Determines the consensus sequence from the motifs.
- `Score(Motifs)`: Calculates the score of the motifs based on the consensus sequence.
- `Pr(kmer, profile)`: Computes the probability of a k-mer given a profile matrix.
- `ProfileMostProbableKmer(text, k, profile)`: Finds the most probable k-mer in a text given a profile matrix.
- `Normalize(Probabilities)`: Normalizes a dictionary of probabilities so that they sum to 1.
- `WeightedDie(Probabilities)`: Simulates a weighted die roll based on a dictionary of probabilities.
- `ProfileGeneratedKmer(Text, profile, k)`: Generates a k-mer from text based on a profile matrix.
- `RepeatedGibbsSampler(Dna, k, t, N, times)`: Repeatedly runs the Gibbs Sampler algorithm to find the best motifs.
- -------------------------------------------------------------------------------------------------------------------------
example usage:
from gibbs_sampler import RepeatedGibbsSampler
Sample DNA sequences

sequences = [
"TCGGGGGTTTTT",
"CCGGTGACTTAC",
"ACGGGGATTTTC",
"TTGGGGACTTTT",
"AAGGGGACTTCC"
]
Find 6bp motifs with 100 iterations and 20 restarts

best_motifs = RepeatedGibbsSampler(
Dna=sequences,
k=6,
t=5,
N=100,
times=20
)

print("Best motifs found:", best_motifs)
print("Consensus sequence:", Consensus(best_motifs))

this shows:
    1. Basic import of the main function
    2. Example DNA input (5 sequences)
    3. Typical parameter values
    4. Simple output of results

-----------------------------------------------------------------------------------------------------------------------------
references:

    Lawrence et al. (1993) Science 262(5131):208-214

    Bioinformatics Algorithms (Compeau & Pevzner, 2015)

    UC San Diego CSE 181 Course Materials
