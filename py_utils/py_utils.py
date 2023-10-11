# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 21:27:56 2019

@author: asus
"""
import matplotlib.pyplot as plt
import seaborn as sns
import tensorflow as tf
import numpy as np
import random
import pandas as pd
import psutil
from sklearn import metrics
from time import strftime


def memory():
    info = psutil.virtual_memory()
    print("可用内存：{} MB".format(info.free/1024/1024))


def im2col(input_data, filter_h, filter_w, stride=1, pad=0):
    """

    Parameters
    ----------
    input_data : 由(数据量, 通道, 高, 长)的4维数组构成的输入数据
    filter_h : 卷积核的高
    filter_w : 卷积核的长
    stride : 步幅
    pad : 填充

    Returns
    -------
    col : 2维数组
    """
    N, C, H, W = input_data.shape
    out_h = (H + 2*pad - filter_h)//stride + 1
    out_w = (W + 2*pad - filter_w)//stride + 1

    img = np.pad(input_data, [(0,0), (0,0), (pad, pad), (pad, pad)], 'constant')
    col = np.zeros((N, C, filter_h, filter_w, out_h, out_w))

    for y in range(filter_h):
        y_max = y + stride*out_h
        for x in range(filter_w):
            x_max = x + stride*out_w
            col[:, :, y, x, :, :] = img[:, :, y:y_max:stride, x:x_max:stride]

    col = col.transpose(0, 4, 5, 1, 2, 3).reshape(N*out_h*out_w, -1)
    return col

def convolution(X, weights, bias=0, stride=1, pad=0):
    '''
        功能：计算卷积
        X: NHWC顺序
        weights: HWCN顺序
    '''
    X = X.transpose(0,3,1,2) # 改成NCHW
    weights = weights.transpose(3,2,0,1) # 改成NCHW
    FN, C, FH, FW = weights.shape
    N, C, H, W = X.shape

    out_h = int(1 + (H + 2*pad - FH) / stride)
    out_w = int(1 + (W + 2*pad - FW) / stride)

    col = im2col(X, FH, FW, stride, pad)
    col_W = weights.reshape(FN, -1).T # 卷积核的展开
    out = np.dot(col, col_W) + bias

    out = out.reshape(N, out_h, out_w, -1)

    return out

def onehot2seq(onehot):
    maps = {
        '[1. 0. 0. 0.]': 'A',
        '[0. 1. 0. 0.]': 'T',
        '[0. 0. 1. 0.]': 'G',
        '[0. 0. 0. 1.]': 'C',
        '[0. 0. 0. 0.]': 'N'
    }
    seq = ''
    for one in onehot:
        seq += maps[str(one)]
    return seq


#########################################
import colorsys
import random
 
def get_n_hls_colors(num):
    hls_colors = []
    i = 0
    step = 360.0 / num
    while i < 360:
        h = i
        s = 90 + random.random() * 10
        l = 50 + random.random() * 10
        _hlsc = [h / 360.0, l / 100.0, s / 100.0]
        hls_colors.append(_hlsc)
        i += step
 
    return hls_colors
 
def ncolors(num):
    rgb_colors = []
    if num < 1:
        return rgb_colors
    hls_colors = get_n_hls_colors(num)
    for hlsc in hls_colors:
        _r, _g, _b = colorsys.hls_to_rgb(hlsc[0], hlsc[1], hlsc[2])
        r, g, b = [int(x * 255.0) for x in (_r, _g, _b)]
        rgb_colors.append([r, g, b])
 
    return rgb_colors

def color(value):
    digit = list(map(str, range(10))) + list("ABCDEF")
    if isinstance(value, tuple):
        string = '#'
        for i in value:
            a1 = i // 16
            a2 = i % 16
            string += digit[a1] + digit[a2]
        return string
    elif isinstance(value, str):
        a1 = digit.index(value[1]) * 16 + digit.index(value[2])
        a2 = digit.index(value[3]) * 16 + digit.index(value[4])
        a3 = digit.index(value[5]) * 16 + digit.index(value[6])
        return (a1, a2, a3)

def get_random_color(num):
    '''
        生成num个十六进制颜色代码
    '''
    return list(map(lambda x: color(tuple(x)), ncolors(num)))


#########################################

def subset(alist, idxs):
    '''
        用法：根据下标idxs取出列表alist的子集
        alist: list
        idxs: list
    '''
    sub_list = []
    for idx in idxs:
        sub_list.append(alist[idx])

    return sub_list


#########################################

def split_list(alist, group_num=4, shuffle=True, retain_left=False):
    '''
        用法：将alist切分成group个子列表，每个子列表里面有len(alist)//group个元素
        shuffle: 表示是否要随机切分列表，默认为True
        retain_left: 若将列表alist分成group_num个子列表后还要剩余，是否将剩余的元素单独作为一组
    '''

    _index = list(range(len(alist))) # 保留下标

    # 是否打乱列表
    if shuffle: 
        index = _index.copy()
        random.shuffle(index) # shuffle会改变index的内容，因此需要设置保留下标_index
    
    elem_num = len(alist) // group_num # 每一个子列表所含有的元素数量
    sub_lists = {}
    
    # 取出每一个子列表所包含的元素，存入字典中
    for idx in range(group_num):
        start, end = idx*elem_num, (idx+1)*elem_num
        sub_lists['set'+str(idx)] = subset(alist, index[start:end])
    
    # 是否将最后剩余的元素作为单独的一组
    if retain_left and group_num * elem_num != len(index): # 列表元素数量未能整除子列表数，需要将最后那一部分元素单独作为新的列表
        sub_lists['set'+str(idx+1)] = subset(alist, index[end:])
    
    return sub_lists


#########################################

def form_to_dict(data):
    data_list = data.split('\n')
    data_dict = {}
    for e in data_list:
        s = e.split(':')
        data_dict[s[0]] = s[1].lstrip()
    return data_dict

def ctime():
    print("\n---Time: ",strftime('%Y/%m/%d %H:%M:%S'))


### 绘制热图
def heatmap(data,xLabel,yLabel,title,result_path):
    f, ax = plt.subplots(figsize=(8,5))
    ax = sns.heatmap(data,cmap = 'YlGnBu',ax=ax,annot=True,fmt ='0.3f',
                     xticklabels=xLabel,yticklabels=yLabel)
    
    plt.xlabel('conv2_filter',fontsize=20)
    plt.ylabel('conv1_filter',fontsize=20)
    plt.title(title,fontsize=30)
    plt.savefig(result_path+title+'.pdf')
    

### 读取MNIST数据集
def mnist_data():
    mnist = tf.keras.datasets.mnist
    
    (x_train, y_train), (x_test, y_test) = mnist.load_data()
    x_train, x_test = x_train / 255.0, x_test / 255.0
    
    # add a channels dimension
    x_train = x_train[..., tf.newaxis]
    x_test = x_test[..., tf.newaxis]
    return x_train,y_train,x_test,y_test


def cal_costs(cost_str):
    total_cost = 0
    lines = []
    for l in cost_str.split('\n'):
        if "：" in l:
            line = l.split('：')
            line[1] = float(line[1])
            lines.append(line)

    # print(lines)

    lines = pd.DataFrame(lines, columns=['类别', '金额'])

    total = lines.groupby('类别').sum()

    total.loc['吃'] = total.loc['食堂'] + total.loc['额外吃']
    total.loc['其他'] = lines['金额'].sum() - total.loc['吃']
    
    print(total)
    print("总花销：",lines['金额'].sum())
    return lines, total
