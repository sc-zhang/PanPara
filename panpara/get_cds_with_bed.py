#!/usr/bin/env python
import argparse
import os


def get_opts():
    groups = argparse.ArgumentParser()
    groups.add_argument('-b', '--bed', help='bed file with 4 columns, \"Chr\tstart_pos\tend_pos\tgene_name\"',
                        required=True)
    groups.add_argument('-c', '--cds', help='directory of all cds files, all cds files need end with \".cds\"',
                        required=True)
    groups.add_argument('-o', '--output', help='output cds', required=True)
    return groups.parse_args()


def get_seq_with_list(in_dir, in_list, out_fa):
    fa_db = {}
    for in_fa in os.listdir(in_dir):
        if not in_fa.endswith('.cds'):
            continue
        in_fa = os.path.join(in_dir, in_fa)
        with open(in_fa, 'r') as fin:
            for line in fin:
                if line[0] == '>':
                    id = line.strip()[1:]
                    fa_db[id] = []
                else:
                    fa_db[id].append(line.strip())

    for id in fa_db:
        fa_db[id] = ''.join(fa_db[id])

    with open(in_list, 'r') as fin:
        with open(out_fa, 'w') as fout:
            for line in fin:
                if line[0] == '#':
                    continue
                else:
                    tig = line.strip().split()[3]
                    if tig not in fa_db:
                        continue
                    fout.write(">%s\n%s\n" % (tig, fa_db[tig]))


if __name__ == "__main__":
    opts = get_opts()
    in_list = opts.bed
    in_dir = opts.cds
    out_fa = opts.output
    get_seq_with_list(in_dir, in_list, out_fa)
