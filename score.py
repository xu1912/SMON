#!/usr/bin/python

##Calculate the combined score.

import sys
import os
from score_cal import compute_combined_score_protein_protein, compute_combined_score_protein_protein_new, compute_combined_score_protein_protein_new_full

of=sys.argv[1]
rs=open(of, 'r')
lines=rs.readlines()
ouf=open(of.split(".")[0]+"_a.txt", 'w')
n=0
for line in lines:
        data=line.strip().split("\t")
        ouf.write(line.strip())
        if n==0:
                n += 1
                ouf.write('\tcal_old\tcal_new\n')
                continue

        cor_a = abs(float(data[2]))
        neighborhood = float(data[3])
        neighborhood_transferred = 0
        fusion = float(data[4])
        cooccurence = float(data[5])
        homology = float(data[6])
        coexpression = float(data[7])
        coexpression_transferred = 0
        experiments = float(data[8])
        experiments_transferred = 0
        database = float(data[9])
        database_transferred = 0
        textmining = float(data[10])
        textmining_transferred = 0

        #apply prior correction
        combined_score = compute_combined_score_protein_protein(neighborhood, neighborhood_transferred, fusion, cooccurence, homology, coexpression, coexpression_transferred, experiments, experiments_transferred, database, database_transferred, textmining, textmining_transferred)

        combined_score_new = compute_combined_score_protein_protein_new_full(cor_a,neighborhood, fusion, cooccurence, homology, coexpression, experiments, database, textmining)

        ouf.write('\t'+str(combined_score)+'\t'+str(combined_score_new)+'\n')


rs.close()
ouf.close()
