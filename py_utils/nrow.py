#!/usr/bin/env python

import os
import sys
import argparse
import gzip
import time

sys.path.append('/public/home/fan_lab/wangjie/utils/py_utils/nrow_exe')

from nrow_exe import get_nrow

def timer(func):
    '''
        Run function and report used time
    '''
    def wrapper():
        start_time = time.time()

        func()

        end_time = time.time()
        used_time = end_time - start_time
        print('Total time: %s seconds'%(used_time))
    return wrapper

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """
        可以将参数的默认选项打印在帮助文档中
    """
    pass

def check_files(*args):
    for afile in args:
        if not os.path.exists(afile):
            raise Exception("Error! File `%s` not exists! Please check it!"%(afile))

def get_args():
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument('-i', "--input-file", help="Input file to count number of rows", required=True)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    return args

@timer
def main():
    args = get_args()
    nrow = get_nrow(args.input_file)
    print('File %s has %s rows'%(args.input_file, nrow))

if __name__ == "__main__":
    main()
