#!/usr/bin/env python

import argparse
import numpy as np
import os
import pandas as pd
import scipy.stats as s
import subprocess
import sys
from capy import seq

def parse_args():
	# parse args
	parser = argparse.ArgumentParser(description = "Get het site coverage from MuTect 1 callstats file")
	parser.add_argument("-c", required = True, help = "Path to callstats file", metavar = "callstats_in")
	parser.add_argument("-s", required = True, help = "Path to GATK-formatted SNP site file", metavar = "snplist_in")
	parser.add_argument("-r", required = True, help = "Path to reference FASTA (directory must contain FASTA index)", metavar = "ref_in")
	parser.add_argument("-o", required = True, help = "Het coverage file prefix ('tumor'/'normal' appended)", metavar = "output_prefix")
	parser.add_argument("-g", help = "Output genotype file", action = "store_true")
	parser.add_argument("--af_lb", help = "Lower bound on beta distribution AF interval", default = 0.4, type = float, metavar = "lowerbound")
	parser.add_argument("--af_ub", help = "Upper bound on beta distribution AF interval", default = 0.6, type = float, metavar = "upperbound")
	parser.add_argument("--dens", help = "Beta distribution density threshold to consider a site heterozygous in the normal.", default = 0.7, type = float, metavar = "cutoff")

	args = parser.parse_args()

	# check args
	if not 0 < args.dens < 1:
		raise ValueError("Density threshold must be between 0 and 1!")

	if not 0 < args.af_lb < 1:
		raise ValueError("AF lower bound must be between 0 and 1!")
	if not 0 < args.af_ub < 1:
		raise ValueError("AF upper bound must be between 0 and 1!")
	if args.af_ub <= args.af_lb:
		raise ValueError("AF lower bound must be less than or equal to upper bound!")

	if not os.path.exists(args.c):
		raise FileNotFoundError("Callstats file not found!")
	if not os.path.exists(args.s):
		raise FileNotFoundError("SNP site file not found!")
	if not os.path.exists(args.r):
		raise FileNotFoundError("Reference fasta file not found!")
	if not os.path.exists(args.r + '.fai'):
		raise FileNotFoundError("Reference fasta index file not found! (Must be <reference.fa>.fai)")

	return args

def hash_altref(DF):
	return (DF.replace(dict(zip(list("ACGT"), range(0, 4))))@[4, 1]).astype(np.uint8)


if __name__ == "__main__":
	args = parse_args()

	# trim callstats (faster to do this on the shell)
	callstats_trimmed = subprocess.Popen("sed '1,2d' {} | cut -f1,2,4,5,26,27,38,39".format(args.c), shell = True, stdout = subprocess.PIPE)

	# load in callstats
	print("Loading callstats file ...", file = sys.stderr)
	CS = pd.read_csv(callstats_trimmed.stdout, sep = "\t",
	  names = ["chr", "pos", "ref", "alt", "t_refcount", "t_altcount", "n_refcount", "n_altcount"],
	  dtype = { "chr" : str, "pos" : np.uint32, "t_refcount" : np.uint32, "t_altcount" : np.uint32, "n_refcount" : np.uint32, "n_altcount" : np.uint32 }
	)
	contig_list = pd.read_csv(args.r + '.fai', sep='\t', usecols = [0], names=["contig"])["contig"].tolist()
	CS["chr"] = CS["chr"].apply(lambda x: contig_list.index(x) + 1).astype(np.uint8)
	CS["gpos"] = seq.chrpos2gpos(CS["chr"], CS["pos"], ref = args.r)
	CS["allele"] = hash_altref(CS.loc[:, ["alt", "ref"]])
	CS = CS.drop(columns = ["alt", "ref"])
	print("{} sites loaded.".format(CS.shape[0]), file = sys.stderr)

	# load in SNP list
	print("Loading SNP list ...", file = sys.stderr)
	H = pd.read_csv(args.s, sep = "\t", comment = "@",
	  names = ["chr", "pos", "x", "y", "allele"],
	  dtype = { "chr" : str, "pos" : np.uint32, "x" : np.uint32, "y" : str, "allele" : str },
	).drop(columns = ["x", "y"])
	H["chr"] = H["chr"].apply(lambda x: contig_list.index(x) + 1).astype(np.uint8)
	H["gpos"] = seq.chrpos2gpos(H["chr"], H["pos"], ref = args.r) 
	H["allele"] = hash_altref(H["allele"].str.extract(r"(.)/(.)"))
	print("{} sites loaded.".format(H.shape[0]), file = sys.stderr)

	# merge
	print("Pulling down SNP site coverage from callstats ...", file = sys.stderr)
	H = H.merge(CS.drop(columns = ["chr", "pos"]), left_on = ["gpos", "allele"], right_on = ["gpos", "allele"], how = "inner")
	print("{} covered SNP sites identified.".format(H.shape[0]), file = sys.stderr)

	# compute which sites in the SNP list are confidently heterozygous in the normal
	# bdens = \int_{af_lb}^{af_ub} df beta(f | n_alt + 1, n_ref + 1)
	H["bdens"] = s.beta.cdf(args.af_ub, H["n_altcount"].values[:, None] + 1, H["n_refcount"].values[:, None] + 1) - \
	s.beta.cdf(args.af_lb, H["n_altcount"].values[:, None] + 1, H["n_refcount"].values[:, None] + 1)

	# compute which sites are confidently homozygous alt. in the normal
	H["bdens_hom"] = 1 - s.beta.cdf(0.95, H["n_altcount"].values[:, None] + 1, H["n_refcount"].values[:, None] + 1)

	# save tumor het coverage at good sites to file (GATK GetHetCoverage format)
	good_idx = H["bdens"] > args.dens
	print("Identified {} high quality het sites in normal.".format(good_idx.sum()), file = sys.stderr)
	H.loc[good_idx, ["chr", "pos", "t_refcount", "t_altcount"]].rename(columns = { "chr" : "CONTIG", "pos" : "POSITION", "t_refcount" : "REF_COUNT", "t_altcount" : "ALT_COUNT" }).to_csv(args.o + ".tumor.tsv", sep = "\t", index = False)

	# save normal het coverage at good sites to file
	H.loc[good_idx, ["chr", "pos", "n_refcount", "n_altcount"]].rename(columns = { "chr" : "CONTIG", "pos" : "POSITION", "n_refcount" : "REF_COUNT", "n_altcount" : "ALT_COUNT" }).to_csv(args.o + ".normal.tsv", sep = "\t", index = False)

	# if requested, save genotype file as TSV (23andme style) 
	if args.g:
		het_idx = good_idx
		hom_idx = H["bdens_hom"] > args.dens
		gen_idx = het_idx | hom_idx
		G = H.loc[gen_idx, ["chr", "pos", "allele"]]

		# add genotype info
		alt_ref = np.array(["A", "C", "G", "T"])[np.c_[(G["allele"].values & 0xC) >> 2, G["allele"].values & 3]]
		alt_ref[hom_idx[gen_idx], 1] = alt_ref[hom_idx[gen_idx], 0] 
		G["genotype"] = np.char.add(alt_ref[:, -1], alt_ref[:, 0])

		# save
		G.drop(columns = ["allele"]).to_csv(args.o + ".genotype.tsv", sep = "\t", index = False)
