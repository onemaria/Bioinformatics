#the more conserved the column, the smaller its entropy. 
# Thus, entropy offers an improved method of scoring motif matrices: 
# the entropy of a motif matrix is defined as the sum of the entropies of its columns

import math

def Entropy(motifs):
    entropy = 0
    for i in range(len(motifs)):
        for j in motifs[i]:
            if j == 0:
                entropy += 0
            else:
                entropy += j*(math.log(j,2))
    return entropy *-1


profile =[
[0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
[0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
[0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
]
print(Entropy(profile))