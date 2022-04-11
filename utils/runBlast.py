#!/usr/bin/env python
import os
import argparse


def get_opts():
    groups = argparse.ArgumentParser()
    groups.add_argument('-q', '--query', help='query fasta', required=True)
    groups.add_argument('-d', '--db', help='reference fasta', required=True)
    groups.add_argument('-p', '--program', help='blast program for use, should be blastn or blastp, default=blastn', default="blastn")
    groups.add_argument('-e', '--evalue', help='evalue, default=1e-3', default="1e-3")
    groups.add_argument('-f', '--format', help='output format of blast, default=6', default="6")
    groups.add_argument('-n', '--num_alignment', help='number of alignment, if is empty, means all', default="")
    groups.add_argument('-o', '--output', help='output blast file', required=True)
    groups.add_argument('-t', '--thread', help='threads for blast, default=6', default="6")
    return groups.parse_args()


def run_blast(qry, ref, prog, evalue, fmt, num_aln, out_blast, threads):
    # Make blast db
    if prog == "blastn":
        db_type = "nucl"
    else:
        db_type = "prot"
    idx_cmd = "makeblastdb -in %s -dbtype %s -out blastdb &> makeblastdb.log"%(ref, db_type)
    if num_aln != "":
        num_aln_cmd = " -num_alignments %s "%num_aln
    else:
        num_aln_cmd = ""
    blast_cmd = "%s -query %s -db blastdb -out %s -evalue %s -outfmt %s %s -num_threads %s &> blast.log"%(prog, qry, out_blast, evalue, fmt, num_aln_cmd, threads)

    print("Running makeblastdb: %s"%idx_cmd)
    ret = os.system(idx_cmd)
    if ret != 0:
        print("Fatal error, makeblastdb failed")
        exit(-1)
    print("Running blast: %s"%blast_cmd)
    ret = os.system(blast_cmd)
    if ret != 0:
        print("Fatal error, blast failed")
        exit(-1)
    print("Finished")


if __name__ == "__main__":
    opts = get_opts()
    qry = opts.query
    ref = opts.db
    prog = opts.program
    evalue = opts.evalue
    fmt = opts.format
    num_aln = opts.num_alignment
    out_blast = opts.output
    threads = opts.thread
    run_blast(qry, ref, prog, evalue, fmt, num_aln, out_blast, threads)
