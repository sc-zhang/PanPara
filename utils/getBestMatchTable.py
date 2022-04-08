#!/usr/bin/env python
import argparse
import os
import sys
from .UnionFind import UnionFind
from .runMCscanX import run_mcscanx


def get_opt():
    group = argparse.ArgumentParser()
    group.add_argument('-1', '--bed1', help="Input first bed", required=True)
    group.add_argument('-2', '--bed2', help="Input second bed", required=True)
    group.add_argument('-a', '--table', help="Paralog table, default=\"\"", default="")
    group.add_argument('-b', '--blast', help="Input blast file", required=True)
    group.add_argument('-l', '--blast2', help="Self blast of query, if bed1 is different with bed2, this parameter is required, default=\"\"", default="")
    group.add_argument('-d', '--iden', help="Identity threshold, default=0.8", default=0.8, type=float)
    group.add_argument('-c', '--coverage', help="The threshold of alignment coverage, default=0.8", default=0.8, type=float)
    group.add_argument('-o', '--output', help="Output directory", required=True)
    return group.parse_args()


def getGeneNo(gene_db, idx_db, len_db, bedfile):
    bed_list = []
    with open(bedfile, 'r') as fin:
        for line in fin:
            data = line.strip().split()
            bed_list.append([data[0], int(data[1]), int(data[2]), data[3]])
    idx = 0
    for chrn, sp, ep, gn in sorted(bed_list):
        # remove duplication genes
        if gn in gene_db:
            continue
        gene_db[gn] = [idx, chrn, sp, ep]
        idx_db[idx] = gn
        len_db[idx] = abs(ep-sp)+1
        idx += 1


def get_best_match_table(bed1, bed2, tbl, blast, blast2, iden_threshold, cov_threshold, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    if bed1 == bed2:
        is_self = True
    else:
        is_self = False
    
    # Self compare is only used for construct first reference
    if is_self and tbl != "":
        print("Fatal error, cannot append infomations to self compare")
        sys.exit(-1)
    # Sorted gene with position, and generate index dict, length dict
    print("Loading bed1")
    bed_db1 = {}
    idx_db1 = {}
    len_db1 = {}
    getGeneNo(bed_db1, idx_db1, len_db1, bed1)
    
    if not is_self:
        print("Loading bed2")
        bed_db2 = {}
        idx_db2 = {}
        len_db2 = {}
        getGeneNo(bed_db2, idx_db2, len_db2, bed2)
    
    # Get score between two genes, and fill the matrix
    match_db = {}
    match_score = {}
    col_dir = outdir.split('.')
    col_dir[-1] = "msx"
    col_out_dir = '.'.join(col_dir)
    col_score_db = run_mcscanx(bed1, bed2, blast, blast2, col_out_dir)
    
    print("Scanning blast")
    with open(blast, 'r') as fin:
        for line in fin:
            data = line.strip().split()
            iden = float(data[2])
            bs = float(data[-1])
            al = int(data[3])
            qn = data[0]
            rn = data[1]
            c_score = 1
            if qn in col_score_db:
                c_score *= col_score_db[qn]
            if rn in col_score_db:
                c_score *= col_score_db[rn]
            
            score = bs*1.0/al*c_score

            if is_self:
                # Skip gene match itself and identity lower than threshold and genes from other samples
                if qn == rn or iden < iden_threshold or qn not in bed_db1 or rn not in bed_db1:
                    continue
                fi = bed_db1[qn][0]
                si = bed_db1[rn][0]

                # If the alignment coverage the gene lower than coverage threshold, skip
                if al*2.0/(len_db1[si]+len_db1[fi]) < cov_threshold:
                    continue

                # Get the highest score and fill in the matrix
                if si not in match_score or score > match_score[si]:
                    match_score[si] = score
                    # Get best matches
                    match_db[si] = set()
                    match_db[si].add(fi)
                elif score == match_score[si]:
                    match_db[si].add(fi)
                
            else:
                # Skip gene pairs with lower identity than threshold
                if iden < iden_threshold:
                    continue
                
                # Skip gene pairs from other samples
                if (qn not in bed_db1 and qn not in bed_db2) or (rn not in bed_db1 and rn not in bed_db2):
                    continue
                
                # Skip gene pairs from same sample
                if (qn in bed_db1 and rn in bed_db1) or (qn in bed_db2 and rn in bed_db2):
                    continue
                
                if qn in bed_db1:
                    fi = bed_db1[qn][0]
                    si = bed_db2[rn][0]
                else:
                    si = bed_db2[qn][0]
                    fi = bed_db1[rn][0]

                # Get the highest score and fill in the matrix
                if si not in match_score or score > match_score[si]:
                    match_score[si] = score
                    match_db[si] = set()
                    match_db[si].add(fi)
                elif score == match_score[si]:
                    match_db[si].add(fi)
    
    # First, scan all queries, and get the best references they matched 
    print("Getting best matches")
    ref_db = {}
    for qry in match_db:
        # For each query, add it to the dict with the reference it best matched
        for idx in match_db[qry]:
            if idx not in ref_db:
                ref_db[idx] = []
            ref_db[idx].append([qry, match_score[qry]])
    
    print("Writing results")
    outfile = os.path.join(outdir, "para.csv")
    new_ref = os.path.join(outdir, "ref.bed")
    bed2_name = bed2.split('/')[-1].split('.')[0]
    
    # If tbl is null, means here is the first time to construct paralog table
    if tbl == "" and (not is_self):
        print("Fatal error, if not self comparison, must afford table constructed by self comparison")
        sys.exit(-1)
    if is_self:
        with open(outfile, 'w') as fout:
            with open(new_ref, 'w') as fref:
                fout.write("#REF,%s\n"%bed2_name)

                # For self comparison, use unionfind to connect paralog table
                ref_cnt = len(bed_db1)
                uf = UnionFind(ref_cnt)
                for ref in ref_db:
                    for qry, _ in ref_db[ref]:
                        uf.union(ref, qry)

                # Connect paralogs                
                para_db = {}
                for ref in idx_db1:
                    gid = uf.find(ref)
                    if gid not in para_db:
                        para_db[gid] = []
                    para_db[gid].append([ref, len_db1[ref]])
                
                para_db_new = {}
                for gid in para_db:
                    tmp = []
                    for ref, _ in sorted(para_db[gid], key = lambda x: x[1], reverse=True):
                        tmp.append(ref)
                    ref = tmp[0]
                    if ref not in para_db_new:
                        para_db_new[ref] = []
                    if len(tmp) > 1:
                        para_db_new[ref] = tmp[1:]

                for ref in sorted(para_db_new):
                    ref_gn = idx_db1[ref]
                    
                    # Write new reference bed
                    fref.write("%s\t%d\t%d\t%s\n"%(bed_db1[ref_gn][1], bed_db1[ref_gn][2], bed_db1[ref_gn][3], ref_gn))
                    # Add paralogs
                    info = [ref_gn]
                    for qry in para_db_new[ref]:
                        qry_gn = idx_db1[qry]
                        info.append(qry_gn)
                    fout.write("%s,"%ref_gn)
                    fout.write("%s\n"%('|'.join(map(str, info))))
    else:
        exists_gn = set()
        with open(outfile, 'w') as fout:
            with open(new_ref, 'w') as fref:
                with open(tbl, 'r') as fin:
                    for line in fin:
                        # Write header
                        if line[0] == '#':
                            col_cnt = len(line.strip().split(','))
                            fout.write("%s,%s\n"%(line.strip(), bed2_name))
                        else:
                            data = line.strip().split(',')
                            ref_gn = data[0]
                            ref = bed_db1[ref_gn][0]
                            # Write new reference bed
                            fref.write("%s\t%d\t%d\t%s\n"%(bed_db1[ref_gn][1], bed_db1[ref_gn][2], bed_db1[ref_gn][3], ref_gn))
                            info = []
                            if ref in ref_db:
                                for qry, _ in sorted(ref_db[ref], key=lambda x:x[1], reverse=True):
                                    qry_gn =idx_db2[qry]
                                    if qry_gn in exists_gn:
                                        continue
                                    exists_gn.add(qry_gn)
                                    info.append(qry_gn)
                            data.append('|'.join(info))
                            fout.write("%s\n"%(','.join(data)))
                
                # The genes not matched with reference should construct paralogs and add to refernce
                with open(blast2, 'r') as fin:
                    nomatch_qry = set()
                    for qry_gn in sorted(bed_db2):
                        qry = bed_db2[qry_gn][0]
                        if qry not in match_db and qry_gn not in exists_gn:
                            nomatch_qry.add(qry_gn)
                       
                    # Get retain queries self comparison
                    col_dir[-1] = "msx_self"
                    col_out_dir = '.'.join(col_dir)
                    col_score_db = run_mcscanx(bed2, bed2, blast2, blast2, col_out_dir)
                    match_score = {}
                    match_db = {}
                    for line in fin:
                        data = line.strip().split()
                        iden = float(data[2])
                        bs = float(data[-1])
                        al = int(data[3])
                        qn = data[0]
                        rn = data[1]
                        c_score = 1
                        if qn in col_score_db:
                            c_score *= col_score_db[qn]
                        if rn in col_score_db:
                            c_score *= col_score_db[rn]

                        score = bs*1.0/al*c_score
                        if qn == rn or iden < iden_threshold or qn not in nomatch_qry or rn not in nomatch_qry:
                            continue
                        fi = bed_db2[qn][0]
                        si = bed_db2[rn][0]

                        if al*2.0/(len_db2[si]+len_db2[fi]) < cov_threshold:
                            continue
                        if si not in match_score or score > match_score[si]:
                            match_score[si] = score
                            match_db[si] = set()
                            match_db[si].add(fi)
                        elif score == match_score[si]:
                            match_db[si].add(fi)
                    
                    # Get best matches
                    qry_db = {}
                    for qry in match_db:
                        for idx in match_db[qry]:
                            if idx not in qry_db:
                                qry_db[idx] = []
                            qry_db[idx].append([qry, match_score[qry]])

                    # For self comparison, use unionfind to connect paralog table
                    qry_cnt = len(nomatch_qry)
                    nomatch_idx = {}
                    nomatch_db = {}
                    idx = 0
                    for qry_gn in sorted(nomatch_qry):
                        qry = bed_db2[qry_gn][0]
                        nomatch_idx[idx] = qry
                        nomatch_db[qry] = idx
                        idx += 1
                    
                    uf = UnionFind(qry_cnt)
                    for qry in qry_db:
                        for next_qry, _ in qry_db[qry]:
                            uf.union(nomatch_db[qry], nomatch_db[next_qry])

                    # Connect paralogs                
                    para_db = {}
                    for idx in range(qry_cnt):
                        gid = uf.find(idx)
                        if gid not in para_db:
                            para_db[gid] = []
                        qry = nomatch_idx[idx]
                        para_db[gid].append([qry, len_db2[qry]])
                    
                    para_db_new = {}
                    for gid in para_db:
                        tmp = []
                        for qry, _ in sorted(para_db[gid], key = lambda x: x[1], reverse=True):
                            tmp.append(qry)
                        qry = tmp[0]
                        para_db_new[qry] = []
                        if len(tmp) > 1:
                            para_db_new[qry] = tmp[1:]
                    
                    for qry in sorted(para_db_new):
                        qry_gn = idx_db2[qry]
                        fref.write("%s\t%d\t%d\t%s\n"%(bed_db2[qry_gn][1], bed_db2[qry_gn][2], bed_db2[qry_gn][3], qry_gn))
                        info = ["" for _ in range(col_cnt)]
                        info[0] = qry_gn
                        tmp = [qry_gn]
                        for next_qry in para_db_new[qry]:
                            next_qry_gn = idx_db2[next_qry]
                            tmp.append(next_qry_gn)
                        info.append('|'.join(tmp))
                        fout.write("%s\n"%(','.join(info)))
    print("Finished")


if __name__ == "__main__":
    opts = get_opt()
    bed1 = opts.bed1
    bed2 = opts.bed2
    tbl = opts.table
    blast = opts.blast
    blast2 = opts.blast2
    iden_threshold = opts.iden
    cov_threshold = opts.coverage
    outdir = opts.output
    get_best_match_table(bed1, bed2, tbl, blast, blast2, iden_threshold, cov_threshold, outdir)
