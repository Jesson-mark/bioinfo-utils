import sys
import string
import json
import itertools as it
from operator import itemgetter
import collections
from argparse import ArgumentParser

### NODE
# copied and modified from /public/home/fan_lab/wangjie/Programs/source_files/mosdepth/scripts/plot-dist.py
# 
# Example
# dist_file=mosdepth/pygmy.1mb.mosdepth.global.dist.txt
# out_tsv=outputs/genome_coverage_dist.tsv
# python scripts/gen_dist_data.py $dist_file -o $out_tsv

def main():
    args = get_args()
    traces = collections.defaultdict(list)
    chroms = collections.OrderedDict()
    chroms["total"] = True

    for f in args.input:
        sample = f.replace(".mosdepth.global.dist.txt", "")
        gen = (x.rstrip().split("\t") for x in open(f))
        for chrom, data in it.groupby(gen, itemgetter(0)):
            if chrom.startswith("GL"):
                continue
            if "Un" in chrom: continue
            if "random" in chrom or "HLA" in chrom: continue
            if chrom.endswith("alt"): continue
            chroms[chrom] = True
            xs, ys = [], []
            v50 = 0
            found = False
            for _, x, y in data:
                y = float(y)
                if y < 0.01:
                    continue
                if not found and y > 0.5:
                    v50 = x
                    found = True
                    print("{}\t{}\t{}\t{:.3f}".format(sample, chrom, x, y))

                xs.append(float(x))
                ys.append(y)

            if len(xs) > 100:
                xs = [x for i, x in enumerate(xs) if ys[i] > 0.02]
                ys = [y for y in ys if y > 0.02]
                if len(xs) > 100:
                    xs = xs[::2]
                    ys = ys[::2]

            traces[chrom].append({
                'xs': [round(x, 3) for x in xs],
                'ys': [round(y, 3) for y in ys],
                'v50' :v50
            })

    try:
        with open(args.output, "w") as outfile:
            outfile.write('chrom\txs\tys\tv50\n')
            for chrom, data in traces.items():
                data = data[0]
                v50 = data['v50']
                for x, y in zip(data['xs'], data['ys']):
                    outfile.write('%s\t%s\t%s\t%s\n'%(
                        chrom, x, y, v50
                    ))
            print('Done')
    except FileNotFoundError:
        sys.exit("ERROR: failed creating output file, does the path exist?")


def get_args():
    parser = ArgumentParser(description="Creates dist data from mosdepth results.")
    parser.add_argument("-o", "--output",
                        default="dist.tsv",
                        help="path and name of output file. Directories must exist.")
    parser.add_argument("input",
                        nargs='+',
                        help="the dist file(s) to use for plotting")
    return parser.parse_args()


if __name__ == '__main__':
    main()

#%%
a = [1,2,3]
b = [4,5,6]

