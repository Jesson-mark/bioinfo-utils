import sys
import gzip
import time
import yaml
import random
import numpy as np
import collections
import subprocess
from scipy.stats import entropy
import logging

def get_logger(logfile = None):
    logger = logging.getLogger()
    logger.setLevel(level=logging.INFO)

    formatter = logging.Formatter('%(levelname)s %(asctime)s %(filename)s %(funcName)s[line:%(lineno)d] : %(message)s')

    # write into file
    if logfile is not None:
        file_handler = logging.FileHandler(logfile, mode='w')
        file_handler.setLevel(level=logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # write to stdout
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    return logger

def run_shell_cmd(cmd):
    with subprocess.Popen(
        cmd, 
        shell=True, 
        bufsize=1, 
        executable = '/usr/bin/bash', 
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, 
        text=True
    ) as proc:
        for line in proc.stdout:
            print(line, end='')

    if proc.returncode != 0:
        sys.exit(proc.returncode)

class ParamParser:
    def __init__(self, params_file):
        self.params_file = params_file
        self._set_params(params_file)
        

    def _set_params(self, params_file):
        with open(params_file, 'r') as ifl:
            params = yaml.load(ifl.read(), Loader=yaml.Loader)

        print(params)

        for k, v in params.items():
            setattr(self, k, v)

def timer(func):
    '''
        Run function and report used time
    '''
    def wrapper():
        start_time = time.time()

        func()

        end_time = time.time()
        used_time = end_time - start_time
        used_mins = used_time / 60
        print('Total time: %.3s minutes'%(used_mins))
    return wrapper

def open_file(afile):
    if afile.endswith('gz'):
        file_reader = gzip.open(afile, 'rt')
    else:
        file_reader = open(afile)
    return file_reader

def get_nrow(afile):
    count = 0
    file_reader = open_file(afile)
    for index, line in enumerate(file_reader):
        count += 1
    return count

def dict_kv_to_str(adict, kv_sep=":", item_sep=";"):
    '''
        kv表示key, value
        kv_sep用于分隔key和value
        item_sep用于分隔不同的item
        example:
            b={'A':2, 'B':4, 'C':5}
            dict_kv_to_str(b) # 'A:2;B:4;C:5'
    '''
    return item_sep.join(['%s%s%s'%(k, kv_sep, v) for k, v in adict.items()])

def cal_gc_content(seq):
    """ Calculate GC Content in a DNA/RNA sequence """
    return round((seq.count('C') + seq.count('G')) / len(seq), 3)

def estimate_shannon_entropy(dna_sequence):
    bases = collections.Counter(dna_sequence)

    probs = [x/sum(bases.values()) for x in bases.values()]
 
    # use scipy to calculate entropy
    entropy_value = entropy(probs, base=2)
 
    return round(entropy_value, 3)

def random_split(alist, n_batch, shuffle=True):
    shuffled_data = alist
    if shuffle:
        random.shuffle(shuffled_data)
    k, m = divmod(len(alist), n_batch)
    batches = (alist[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n_batch))
    return [b for b in batches if b]

def get_base_score(base1, base2, match=0, mismatch=1):
    '''
        若两个碱基相同，则返回match得分，若不同则返回mismtch得分
    '''
    return match if base1 == base2 else mismatch

def cal_edit_distance(seq1, seq2, match=0, mismatch=1, indel=1):
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
            mis_up = score_matrix[j-1,i] + indel
            mis_left = score_matrix[j,i-1] + indel
            score_matrix[j,i] = min(mis_up, mis_left, base_score+score_matrix[j-1,i-1])

    final_score = score_matrix[-1,-1] # 两条序列的最佳匹配得分

    return int(final_score)
