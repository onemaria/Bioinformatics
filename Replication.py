# the Minimum Skew Problem provides us with an approximate location of ori at position 
# in most bacteria nucleotide compositions are asymmetric between the leading strand and the lagging strand: 
# the leading strand contains more guanine (G) and thymine (T), whereas the lagging strand contains more adenine (A) and cytosine (C)

def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions

def ApproximatePatternCount(Pattern, Text, d):
    count = 0
    for i in range(len(Text)-len(Pattern) + 1):
        if HammingDistance(Pattern, Text[i:i+len(Pattern)]) <= d:
            count+=1
    return count

def HammingDistance(p, q):
    hamdistance = 0
    for i in range(len(p)):
        if p[i]!=q[i]:
            hamdistance+=1
    return hamdistance

def MinimumSkew(Genome):
    array = SkewArray(Genome)
    positions = []
    count = 0
    minarray = min(array)
    for i in array:
        if i == minarray:
            positions.append(count)
        count +=1
    return positions

def SkewArray(Genome):
    skew = [0]
    score = {"A":0, "T":0, "C":-1, "G":1}
    for i in range(1,len(Genome)+1):
            skew.append(score[Genome[i-1]] + skew[i-1])
    return skew
