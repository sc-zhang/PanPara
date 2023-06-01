#!/usr/bin/env python
import os
import argparse
from panpara import *
import time


def time_print(info):
    print("\033[32m%s\033[0m %s" % (time.strftime('[%H:%M:%S]', time.localtime(time.time())), info))


def get_opts():
    groups = argparse.ArgumentParser()
    groups.add_argument('-l', '--list', help='list file, each row contain one sample name, first row means reference',
                        required=True)
    groups.add_argument('-s', '--cds',
                        help='cds directory, must exists all samples end with \".cds\", for example: \"sample1.cds\"',
                        required=True)
    groups.add_argument('-b', '--bed',
                        help='bed directory, must exists all samples end with \".bed\" , for example: \"sample1.bed\"',
                        required=True)
    groups.add_argument('-d', '--iden', help="identity threshold, default=0.8", default=0.8, type=float)
    groups.add_argument('-c', '--coverage', help="the threshold of alignment coverage, default=0.8", default=0.8,
                        type=float)
    groups.add_argument('-o', '--output', help="output directory", required=True)
    groups.add_argument('-t', '--thread', help='threads for running some steps, default=6', default="6")
    return groups.parse_args()


def pan_para(in_list, cds_dir, bed_dir, iden, cov, outdir, threads):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    cds_dir = os.path.abspath(cds_dir)
    bed_dir = os.path.abspath(bed_dir)
    ref_bed = ""
    tbl = ""
    ref_cds = ""
    # get sample list
    time_print("Loading samples")
    sample_list = []
    with open(in_list, 'r') as fin:
        for line in fin:
            smp = line.strip()
            if len(smp) != 0:
                sample_list.append(smp)

    for iter in range(1, len(sample_list) + 1):
        # For each iteration, we need enter each directory, that curdir is use for returning current directory
        curdir = os.path.abspath(".")
        time_print("Starting iteration %d" % iter)
        if iter == 1:
            smp = sample_list[0]

            cds1 = os.path.join(cds_dir, "%s.cds" % smp)
            bed1 = os.path.join(bed_dir, "%s.bed" % smp)
            cds1 = os.path.abspath(cds1)
            bed1 = os.path.abspath(bed1)

            iter_path = os.path.join(outdir, "iter1_%s_%s" % (smp, smp))
            iter_path = os.path.abspath(iter_path)
            if not os.path.exists(iter_path):
                os.makedirs(iter_path)
            os.chdir(iter_path)

            blast_file = os.path.join(iter_path, "iter1.blast")
            if not os.path.exists("blast.ok"):
                # Self comparison
                time_print("\trunning blast")
                run_blast(cds1, cds1, 'blastn', '1e-3', '6', '', blast_file, threads)
                os.system("touch blast.ok")
            else:
                time_print("\tblast finished before, skip")

            # Generate base paralog file
            match_dir = os.path.join(iter_path, "match")
            if not os.path.exists("match.ok"):
                time_print("\tgetting paralogs")
                get_best_match_table(bed1, bed1, '', blast_file, '', iden, cov, match_dir)
                os.system("touch match.ok")
            else:
                time_print("\tparalogs get before, skip")
        else:
            smp = sample_list[iter - 1]

            cds2 = os.path.join(cds_dir, "%s.cds" % smp)
            bed2 = os.path.join(bed_dir, "%s.bed" % smp)
            cds2 = os.path.abspath(cds2)
            bed2 = os.path.abspath(bed2)

            iter_path = os.path.join(outdir, "iter%d_ref%d_%s" % (iter, iter - 1, smp))
            iter_path = os.path.abspath(iter_path)
            if not os.path.exists(iter_path):
                os.makedirs(iter_path)
            os.chdir(iter_path)

            blast_file = os.path.join(iter_path, "iter%d.blast" % iter)
            if not os.path.exists("blast.ok"):
                time_print("\trunning blast")
                run_blast(cds2, ref_cds, 'blastn', '1e-3', '6', '', blast_file, threads)
                os.system("touch blast.ok")
            else:
                time_print("\tblast finished before, skip")

            qry_self_blast_file = os.path.join(iter_path, "iter%d_qry_self.blast" % iter)
            if not os.path.exists("qry_self_blast.ok"):
                time_print("\trunning query self blast")
                run_blast(cds2, cds2, 'blastn', '1e-3', '6', '', qry_self_blast_file, threads)
                os.system("touch qry_self_blast.ok")
            # Gerenate paralog file
            match_dir = os.path.join(iter_path, "match")
            if not os.path.exists("match.ok"):
                time_print("\tgetting paralogs")
                get_best_match_table(ref_bed, bed2, tbl, blast_file, qry_self_blast_file, iden, cov, match_dir)
                os.system("touch match.ok")
            else:
                time_print("\tparalogs get before, skip")

        # Get ref.bed, tbl, ref.cds for next iteration
        ref_bed = os.path.join(match_dir, "ref.bed")
        tbl = os.path.join(match_dir, "para.csv")
        ref_cds = os.path.join(match_dir, "ref.cds")

        if not os.path.exists("newcds.ok"):
            # Generate ref.cds for next iteration
            time_print("\twriting new cds")
            get_seq_with_list(cds_dir, ref_bed, ref_cds)
            os.system("touch newcds.ok")
        else:
            time_print("\tnew cds generated before, skip")

        os.chdir(curdir)

    # Get final result
    time_print("Getting final result")
    final_tbl = os.path.join(outdir, "final.csv")
    os.system("cp %s %s" % (tbl, final_tbl))

    time_print("Finished")


if __name__ == '__main__':
    opts = get_opts()
    in_list = opts.list
    cds_dir = opts.cds
    bed_dir = opts.bed
    iden = opts.iden
    cov = opts.coverage
    outdir = opts.output
    threads = opts.thread
    pan_para(in_list, cds_dir, bed_dir, iden, cov, outdir, threads)
