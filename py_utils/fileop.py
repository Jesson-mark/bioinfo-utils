# 这里存放的是与文件操作有关的函数
import os
import numpy as np
import pandas as pd
import glob

func_list = ['write_str',
             'write_1d_list',
             'write_2d_list',
             'read_lines',
             'read_fasta',
             'write_fasta'
             ]

def pmkdir(path):
    '''
        Make path and return True if it does not exist. 
        If this path exists, return False
    '''
    if not os.path.exists(path):
        os.makedirs(path)
        return True
    else:
        return False

def read_str(filename):
    with open(filename, 'r') as ifl:
        strs = ifl.read()
    
    return strs

def write_str(data,filename):
    with open(filename,'a') as ofl:
        ofl.write(data)

def print_1d_list(data):
    for l in data:
        print(l)

def write_1d_list(data,filename, sep='\t'):
    '''
        将1维列表写入文件中，以sep分隔，即写入到一行中
        sep="\t"
    '''
    with open(filename,'a') as ofl:
        for sub in data:
            if type(sub) == str:
                ofl.write(sub+sep)
            else:
                ofl.write(str(np.round(sub,3))+sep)
        ofl.write('\n')

def write_2d_list(data, filename, sep='\t', round=False):
    ''' 将2维列表按行写入文件，列表中元素类型为数字或字符串 '''
    with open(filename,'a') as ofl:
        for sub in data:
            for ssub in sub:
                if type(ssub) == str:
                    ofl.write(ssub+sep)
                else:
                    if round:
                        ofl.write(str(np.round(ssub,3))+sep)
                    else:
                        ofl.write(str(ssub)+sep)
            ofl.write('\n')


def read_lines(filename, delimiter = '\t', skip_lines = 0):
    '''
        按行读入文件，返回一个二维列表
    '''
    lines = []
    with open(filename) as ifl:
        for i in range(skip_lines):
            ifl.readline() # 跳过文件的前skip_lines行内容
        for l in ifl.readlines():
            ll = l.strip().split(delimiter)
            lines.append(ll)
    return lines

def read_cols(filename, delimiter = '\t', skip_lines = 0, col_type = str, skip_strs=None, pd_df=True, header=True, label=None):
    '''
        按列读取文件，每一列读成一个子列表
        delimiter：分隔符
        skip_lines：跳过文件的前几行
        col_type：指定每一列数据的数据类型
        pd_df: 指定是否要将结果存为pandas的DataFrame形式
        header： 指定是否将文件首行作为DataFrame的头
    '''
    cols = []
    
    # 先确定文件有几列
    with open(filename) as ifl:
        ncol = len(ifl.readline().strip().split(delimiter))
        for i in range(ncol):
            cols.append([])
            
    # 读取文件的所有内容，每一列读成一个列表
    with open(filename) as ifl:
        for i in range(skip_lines):
            ifl.readline() # 跳过文件的前skip_lines行内容

        for l in ifl.readlines():
            if skip_strs is not None: # 跳过以skip_strs开头的行
                if l.startswith(skip_strs):
                    continue
            line = l.strip().split(delimiter)
            ncol = len(line)
            for i in range(ncol):
                cols[i].append(col_type(line[i]))

    # 将文件头加进去
    if header:
        with open(filename) as ifl:    
            header_content = ifl.readline().strip().split(delimiter)
            cols = [header_content] + cols
    
    # 是否返回pd.DataFrame数据类型
    if pd_df:
        new_cols = {}
        for i in range(len(header_content)):
            new_cols[header_content[i]] = cols[i+1]        
        cols = pd.DataFrame(new_cols)

        # 是否将label设置为单独的一列
        if label:
            cols['label'] = label

    return cols

#####读取fasta文件
def read_fasta(filename): 
    ifl=open(filename,'rt')
    iflst=ifl.readlines()
    ifl.close()         

    seqlist=[]
    aseq = []
    titstr = ''
    seqstr=''
    for i in iflst:
        i = i.strip()
        if i[0]=='>':
            titstr = i
            if seqstr!='':
                aseq.append(pretitstr)
                aseq.append(seqstr)
                seqlist.append(aseq)                
                seqstr = ''
                aseq = []
        else:
            seqstr += i
            pretitstr = titstr
    aseq = [titstr, seqstr] 
    seqlist.append(aseq)
    return seqlist

#####将序列与其序列名写入fasta文件，输入为一个二维列表，其内的每一个子元素为[序号，序列]
def write_fasta(data,filename):
    ofl = open(filename,'w')
    for seq in data:
        ofl.writelines([seq[0],'\n',seq[1],'\n'])
    ofl.close()


