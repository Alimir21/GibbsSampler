import random

# Gibbs Sampler
def GibbsSampler(Dna, k, t, N):
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs
    for j in range(N):
        i = random.randint(0, t - 1)
        profile = ProfileWithPseudocounts([Motifs[m] for m in range(t) if m != i])
        Motifs[i] = ProfileGeneratedKmer(Dna[i], profile, k)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

# Generate random motifs
def RandomMotifs(Dna, k, t):
    Motifs = []
    for i in range(t):
        start_index = random.randint(0, len(Dna[i]) - k)
        Motifs.append(Dna[i][start_index:start_index + k])
    return Motifs

# Count with pseudocounts
def CountWithPseudocounts(Motifs):
    count = {nucleotide: [1] * len(Motifs[0]) for nucleotide in 'ACGT'}
    for motif in Motifs:
        for i, nucleotide in enumerate(motif):
            count[nucleotide][i] += 1
    return count

# Profile with pseudocounts
def ProfileWithPseudocounts(Motifs):
    counts = CountWithPseudocounts(Motifs)
    t = len(Motifs) + 4
    profile = {nucleotide: [count / t for count in counts[nucleotide]] for nucleotide in counts}
    return profile

# Consensus sequence
def Consensus(Motifs):
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    consensus = ""
    for i in range(k):
        max_count = 0
        frequent_symbol = ""
        for symbol in "ACGT":
            if count[symbol][i] > max_count:
                max_count = count[symbol][i]
                frequent_symbol = symbol
        consensus += frequent_symbol
    return consensus

# Score motifs
def Score(Motifs):
    consensus = Consensus(Motifs)
    score = sum(motif[i] != consensus[i] for motif in Motifs for i in range(len(motif)))
    return score

# Probability of a k-mer
def Pr(kmer, profile):
    probability = 1.0
    for i, nucleotide in enumerate(kmer):
        probability *= profile[nucleotide][i]
    return probability

# Profile most probable k-mer
def ProfileMostProbableKmer(text, k, profile):
    max_probability = -1
    most_probable_kmer = text[:k]
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        probability = Pr(kmer, profile)
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer
    return most_probable_kmer

# Normalize probabilities
def Normalize(Probabilities):
    total = sum(Probabilities.values())
    return {k: v / total for k, v in Probabilities.items()}

# Weighted die
def WeightedDie(Probabilities):
    p = random.uniform(0, 1)
    cumulative_probability = 0.0
    for kmer, probability in Probabilities.items():
        cumulative_probability += probability
        if p <= cumulative_probability:
            return kmer

# Profile generated k-mer
def ProfileGeneratedKmer(Text, profile, k):
    n = len(Text)
    probabilities = {Text[i:i + k]: Pr(Text[i:i + k], profile) for i in range(n - k + 1)}
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

# Repeated Gibbs Sampler
def RepeatedGibbsSampler(Dna, k, t, N, times):
    BestMotifs = GibbsSampler(Dna, k, t, N)
    for _ in range(times - 1):
        motifs = GibbsSampler(Dna, k, t, N)
        if Score(motifs) < Score(BestMotifs):
            BestMotifs = motifs
    return BestMotifs
