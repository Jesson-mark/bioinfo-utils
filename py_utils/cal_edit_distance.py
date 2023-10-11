#!/usr/bin/env python
#encoding:utf-8

import sys
import numpy as np

def get_base_score(base1, base2, match, mismatch):
    '''
        若两个碱基相同，则返回match得分，若不同则返回mismtch得分
    '''
    return match if base1 == base2 else mismatch

def edit_distance(seq1, seq2, match=0, mismatch=1, indel=1):
    '''
        仅在score_matrix[j,i]这一行与序列比对不同，在这里是取最小值
        可参考该链接：https://users.monash.edu/~lloyd/tildeAlgDS/Dynamic/Edit/#:~:text=The%20simple%20edit%20distance%20algorithm%20would%20normally%20be,a%20distance%20between%20DNA%20sequences%20%28strings%20over%20%7BA%2CC%5D
    '''
    seq1_len = len(seq1)
    seq2_len = len(seq2)
    score_matrix = np.zeros((seq2_len+1, seq1_len+1))

    # 将第一行与第一列置为indel的得分
    for i in range(1, seq1_len+1):
        score_matrix[0, i] = score_matrix[0, i-1] + indel

    for i in range(1, seq2_len+1):
        score_matrix[i, 0] = score_matrix[i-1, 0] + indel

    # 计算score_matrix的其他分值
    for j in range(1, seq2_len+1):
        for i in range(1, seq1_len+1):
            # print(i,j)
            base1 = seq1[i-1]
            base2 = seq2[j-1]
            base_score = get_base_score(base1, base2, match, mismatch)
            mis_up = score_matrix[j-1,i]+indel
            mis_left = score_matrix[j,i-1]+indel
            score_matrix[j,i] = min(mis_up, mis_left, base_score+score_matrix[j-1,i-1])

    final_score = score_matrix[-1,-1] # 两条序列的最佳匹配得分

    return int(final_score)

seq1, seq2 = sys.argv[1:3]
#print('seq1: %s\nseq2: %s'%(seq1,seq2))
print('Edit distance is: %s'%(edit_distance(seq1, seq2)))
