#!/usr/bin/env python

import argparse
from collections import Counter, OrderedDict
import dendropy
from itertools import groupby
import subprocess
import sys
import os
import weblogolib as w
from util_functions import translate, write_to_fasta
partis_path = "./lib/partis" 
sys.path.append(os.path.join(partis_path, 'python'))
import utils as partisutils

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Tabulate the naive sequence posterior probabilities.")
    parser.add_argument(
        'trees_path', type=str,
        help="Path to linearham trees file.")
    parser.add_argument(
        '--partis-yaml-file', type=str,
        help="Path to partis output file for cdr3 annotation.")
    parser.add_argument(
        '--seed-name', type=str,
        help="Name of the seed sequence")
    parser.add_argument(
        '--output-base', type=str, required=True,
        help="The output basename.")

    args = parser.parse_args()

    tree_yielder = dendropy.Tree.yield_from_files(
        files=[args.trees_path],
        schema="newick",
        preserve_underscores=True
    )

    def create_logo(fname):
        with open(fname + ".fasta", "rU") as f: 
            seqs = w.read_seq_data(f)
        data = w.LogoData.from_seqs(seqs)
        subprocess.check_call(("rm " + fname + ".fasta").split())
        
        options = w.LogoOptions()
        options.color_scheme = w.colorscheme.chemistry
        options.unit_name = "probability"
        options.yaxis_label = "Probability"
        options.xaxis_label = "Site Position"
        options.show_fineprint = False
        options.stacks_per_line = 500
        options.tic_length = 10
        
        format = w.LogoFormat(data, options)
        with open(fname + ".png", 'w') as f:
            f.write(w.png_print_formatter(data, format))
    
    naive_seqs = [tree.find_node_with_taxon_label("naive").annotations.get_value("ancestral") for tree in tree_yielder]
    _, annotation_list, cpath = partisutils.read_output(args.partis_yaml_file) 

    if cpath is not None and len(cpath.partitions) > 0:
        annotations = {':'.join(adict['unique_ids']) : adict for adict in annotation_list}  # collect the annotations in a dictionary so they're easier to access
        most_likely_partition = cpath.partitions[cpath.i_best]  # a partition is represented as a list of lists of strings, with each string a sequence id
        sorted_clusters = sorted(most_likely_partition, key=len, reverse=True)
        cluster_annotation = annotations[':'.join(sorted_clusters[0])]
    else:
        cluster_annotation = annotation_list[0]

    aa_cdr3_bounds = (int(cluster_annotation['codon_positions']['v']/3), int((cluster_annotation['codon_positions']['j'] + 3)/3))
    
    aa_naive_seqs = []
    aa_naive_seqs_d, aa_cdr3_seqs_d = {}, {}
    for i, seq in enumerate(naive_seqs):
        aa_naive_seq = translate(seq)
        aa_naive_seqs.append(aa_naive_seq)
        aa_cdr3_seq = aa_naive_seq[aa_cdr3_bounds[0] : aa_cdr3_bounds[1]]
        aa_naive_seqs_d[("naive" + str(i))] = aa_naive_seq
        aa_cdr3_seqs_d[("naive" + str(i))] = aa_cdr3_seq
    full_fname = os.path.join(args.output_base, "aa_naive_seqs_{}".format(args.seed_name))
    cdr3_fname = os.path.join(args.output_base, "cdr3_aa_naive_seqs_{}".format(args.seed_name))
    
    for fname, seqs_d in zip([full_fname, cdr3_fname], [aa_naive_seqs_d, aa_cdr3_seqs_d]):
        write_to_fasta(seqs_d, fname + ".fasta")
        create_logo(fname)

    args.output_base = os.path.join(args.output_base, "aa_naive_seqs")

    aa_naive_seqs_c = Counter(aa_naive_seqs)
    num_trees = len(aa_naive_seqs)
    aa_naive_seqs_d = OrderedDict(
        ("naive_" + str(i) + "_" + str(float(count) / num_trees), seq)
        for i, (seq, count) in enumerate(aa_naive_seqs_c.most_common(None))
    )
    write_to_fasta(aa_naive_seqs_d, args.output_base + ".fasta")

    aa_dna_naive_seqs_d = {}
    for k, g in groupby(naive_seqs, lambda seq: translate(seq)):
        if k in aa_dna_naive_seqs_d:
            aa_dna_naive_seqs_d[k].update(g)
        else:
            aa_dna_naive_seqs_d[k] = Counter(g)

    aa_dna_naive_seqs_map = OrderedDict(
        (k, "\n".join(str(float(count) / num_trees) + "," + dna_seq
                      for dna_seq, count in aa_dna_naive_seqs_d[v].most_common(None)))
        for k, v in aa_naive_seqs_d.iteritems()
    )
    write_to_fasta(aa_dna_naive_seqs_map, args.output_base + ".dnamap")
