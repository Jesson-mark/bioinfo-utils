#!/usr/bin/env python

import os
import sys
import argparse

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    '''
        可以将参数的默认选项打印在帮助文档中
    '''
    pass

def check_files(*args):
    for afile in args:
        if not os.path.exists(afile):
            raise Exception("Error! File `%s` not exists! Please check it!"%(afile))

def get_args():
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("-p", "--project-dir", help="Project dir", required=True)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    return args

def main():
    args = get_args()

if __name__ == "__main__":
    main()

