#!/usr/bin/env python

import gzip

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


