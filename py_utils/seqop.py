
##### 取出DNA序列的反向互补序列
def reverse_complement(seq):
    '''
        seq: string
        使用字符串的translate方法来取出互补序列
    '''
    trans_table = str.maketrans('ACGTacgt', 'TGCAtgca') # 制作映射表
    complement_seq = seq.translate(trans_table) # 取出互补序列
    rev_comp_seq = complement_seq[::-1] # 将互补序列反向
    
    return rev_comp_seq


def reverse_complement_dict(seq):
    '''
        seq: string
        这个函数使用字典映射来取出互补序列
    '''
    comp_dict = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'a': 't',
        't': 'a',
        'g': 'c',
        'c': 'g',
        }
    complement_seq = [comp_dict[base] for base in seq]
    complement_seq = ''.join(complement_seq)
    rev_comp_seq = complement_seq[::-1] # 将互补序列反向
    
    return rev_comp_seq

def shift_sequence(seq, shift_num):
    ''' 
        将序列的两端shift个字符替换为N 
        shift > 0，替换序列右端shift个字符为N
        shift < 0，替换序列左端shift个字符为N
    '''
    shift_range = range(-shift_num, shift_num+1)

    new_seqs = []
    for shift in shift_range:
        if shift > 0:
            new_seq = seq[:(-shift)] + shift * 'N'
        elif shift < 0:
            new_seq = abs(shift) * 'N' + seq[(-shift):]
        else:
            new_seq = seq
        new_seqs.append(new_seq)
    return new_seqs


### 将序列扩充为2000bp
def seq_to_2000(seq): 
    if len(seq) < 2000:
        ln = int((2000-len(seq))/2)
        lr = 2000 - ln - len(seq)
        full_seq = ln*'N'+seq+lr*'N'
    else:
        ln = int((len(seq)-2000)/2)
        full_seq = seq[ln:(2000+ln)]
    return full_seq


def extend_seq(seq, target_length): 
    '''
        将序列扩增或削减为target_length长度，两端补N
    '''
    if len(seq) < target_length:
        ln = int((target_length-len(seq))/2)
        lr = target_length - ln - len(seq)
        full_seq = ln*'N'+seq+lr*'N'
    else:
        ln = int((len(seq)-target_length)/2)
        full_seq = seq[ln:(target_length+ln)]
    return full_seq


#####每条序列转化为三维数组
def str_to_one_hot(string):
    one_hot = []
    switcher = {
        'A': [[1],[0],[0],[0]],
        'a': [[1],[0],[0],[0]],
        
        'T': [[0],[1],[0],[0]],
        't': [[0],[1],[0],[0]],
        
        'C': [[0],[0],[1],[0]],
        'c': [[0],[0],[1],[0]],
        
        'G': [[0],[0],[0],[1]],
        'g': [[0],[0],[0],[1]],
        
        'N': [[0],[0],[0],[0]]
    }
    for i in range(len(string)):
        one_hot.append(switcher.get(string[i]))
    return one_hot



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
def write_fasta(data, filename):
    ofl = open(filename,'w')
    for seq in data:
        ofl.writelines([seq[0],'\n',seq[1],'\n'])
    ofl.close()



#####三维one hot转化为序列
def array_to_seq(array):
    switcher = {
        '[[0]\n [0]\n [0]\n [0]]':'N',
        '[[1]\n [0]\n [0]\n [0]]':'A',
        '[[0]\n [1]\n [0]\n [0]]':'T',
        '[[0]\n [0]\n [1]\n [0]]':'C',        
        '[[0]\n [0]\n [0]\n [1]]':'G'
        }
    k=0
    seq=[]
    for a in array:
        ss=''
        k+=1
        print(k)
        for s in a:
            ss+=switcher.get(str(s))
        seq.append(ss.replace('N',''))
    return seq


