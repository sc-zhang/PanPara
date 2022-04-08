#!/usr/bin/env python
import os
import argparse

#get_best_match_table(bed1, bed2, tbl, blast, blast2, iden_threshold, cov_threshold, outdir)
def get_opts():
    groups = argparse.ArgumentParser()
    groups.add_argument('-1', '--bed1', help="Input first bed", required=True)
    groups.add_argument('-2', '--bed2', help="Input second bed", required=True)
    groups.add_argument('-b', '--blast', help="Input blast file", required=True)
    groups.add_argument('-l', '--blast2', help="Self blast of query, if bed1 is different with bed2, this parameter is required, default=\"\"", default="")
    groups.add_argument('-o', '--output', help='output blast file', required=True)
    groups.add_argument('-t', '--thread', help='threads for blast, default=6', default="6")
    return groups.parse_args()


def load_bed(bed, bed_list):
    with open(bed, 'r') as fin:
        for line in fin:
            data = line.strip().split()
            bed_list.append([data[0][-3:], data[3], data[1], data[2]])


def load_blast(blast, blast_list):
    with open(blast, 'r') as fin:
        for line in fin:
            blast_list.append(line.strip())


def run_mcscanx(bed1, bed2, blast1, blast2, out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    print("Loading bed")
    bed_list = []
    load_bed(bed1, bed_list)
    if bed1 != bed2:
        load_bed(bed2, bed_list)

    print("Loading blast")
    blast_list = []
    load_blast(blast1, blast_list)
    if bed1 != bed2:
        load_blast(blast2, blast_list)

    m_pre = out_dir+"/xyz"
    m_gff = m_pre+'.gff'
    m_blast = m_pre+'.blast'

    print("Writing gff and blast")
    with open(m_gff, 'w') as fout:
        for info in bed_list:
            fout.write("%s\n"%('\t'.join(info)))
    
    with open(m_blast, 'w') as fout:
        fout.write("%s\n"%('\n'.join(blast_list)))
    
    print("Running MCScanX")
    cmd = "MCScanX "+m_pre
    print("\tRunning: %s"%cmd)
    os.system(cmd)

    col_file = m_pre+".collinearity"
    col_db = {}
    with open(col_file, 'r') as fin:
        for line in fin:
            if line[0] == '#':
                if line.startswith("## Alignment"):
                    score = float(line.strip().split()[3].split('=')[1])
            else:
                data = line.strip().split()
                id1 = data[2]
                id2 = data[3]
                if id1 not in col_db or col_db[id1]<score:
                    col_db[id1] = score
                if id2 not in col_db or col_db[id2]<score:
                    col_db[id2] = score

    print("Finished")
    return col_db


if __name__ == "__main__":
    opts = get_opts()
    bed1 = opts.bed1
    bed2 = opts.bed2
    blast1 = opts.blast
    blast2 = opts.blast2
    out_dir = opts.output
    run_mcscanx(bed1, bed2, blast1, blast2, out_dir)
