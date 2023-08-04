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
	parser.add_argument("-r", required = True, help = "Path to reference FASTA (directory must contain FASTA index)", metavar = "ref_in")
	parser.add_argument("-o", required = True, help = "Het coverage file prefix ('tumor'/'normal' appended)", metavar = "output_prefix")
	parser.add_argument("-s", help="Path to GATK-formatted SNP site file", metavar="snplist_in")
	parser.add_argument("-g", help = "Output genotype file", action = "store_true")

	genotyper_parser = parser.add_mutually_exclusive_group(required=False)
	genotyper_parser.add_argument("--use_pod_genotyper", dest="use_pod_genotyper",help = "Use posterior odds method for genotyping", action = "store_true")
	genotyper_parser.add_argument("--use_beta_density", dest="use_pod_genotyper", help = "Use beta distribution density for genotyping", action = "store_false")
	genotyper_parser.add_argument("--use_tonly_genotyper", help = "Genotype in single sample mode", action = "store_true")
	parser.add_argument("--pod_min_depth", type=int, default=10,
						help="Any position with total normal coverage below this threshold will not be considered for genotyping")
	parser.add_argument("--min_tumor_depth", type=int, default=1,
						help="Any position with total tumor coverage below this threshold will be discarded")

	parser.add_argument("--log_pod_threshold", type = float, metavar='threshold', default = 2.5)
	parser.add_argument("--af_lb", help = "Lower bound on beta distribution AF interval", default = 0.4, type = float, metavar = "lowerbound")
	parser.add_argument("--af_ub", help = "Upper bound on beta distribution AF interval", default = 0.6, type = float, metavar = "upperbound")
	parser.add_argument("--dens", help = "Beta distribution density threshold to consider a site heterozygous in the normal.", default = 0.7, type = float, metavar = "cutoff")

	parser.add_argument("--max_frac_mapq0", help = "Any position from callstats with more than this percentage of MAPQ0 reads will be excluded from het coverage analysis.", default = 0.05, type = float, metavar = "mapq0_frac")
	parser.add_argument("--max_frac_prefiltered", help = "Any site with more than this fraction of total reads prefiltered by MuTect will be excluded.", default = 0.10, type = float, metavar = "prefiltered_frac")

	parser.set_defaults(use_pod_genotyper=True)
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

	if not args.log_pod_threshold > 0:
		raise ValueError("Posterior odds threshold must be >0!")

	if not 0 <= args.max_frac_mapq0 <= 1:
		raise ValueError("Max fraction of reads with MAPQ0 at a given position must be between 0 and 1!")
	if not 0 <= args.max_frac_prefiltered <= 1:
		raise ValueError("Max fraction of prefiltered reads must be between 0 and 1!")

	if not os.path.exists(args.c):
		raise FileNotFoundError("Callstats file not found!")
	if args.s is None:
		print("WARNING: without an input SNP site file, all valid het sites will be returned.", file = sys.stderr)
	elif not os.path.exists(args.s):
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
	callstats_trimmed = subprocess.Popen("sed '1,2d' {} | cut -f1,2,4,5,16,17,26,27,38,39".format(args.c), shell = True, stdout = subprocess.PIPE)

	# load in callstats
	print("Loading callstats file ...", file = sys.stderr)
	CS = pd.read_csv(callstats_trimmed.stdout, sep = "\t",
	  names = ["chr", "pos", "ref", "alt", "total_reads", "mapq0_reads", "t_refcount", "t_altcount", "n_refcount", "n_altcount"],
          dtype = { "chr" : str, "pos" : np.uint32, "total_reads" : np.uint32, "mapq0_reads" : np.uint32, "t_refcount" : np.uint32, "t_altcount" : np.uint32, "n_refcount" : np.uint32, "n_altcount" : np.uint32 }
	)
	contig_list = pd.read_csv(args.r + '.fai', sep='\t', usecols = [0], names=["contig"])["contig"].tolist()
	CS["chr"] = CS["chr"].apply(lambda x: contig_list.index(x) + 1).astype(np.uint8)
	CS["gpos"] = seq.chrpos2gpos(CS["chr"], CS["pos"], ref = args.r)
	CS["allele"] = hash_altref(CS.loc[:, ["alt", "ref"]])
	CS = CS.drop(columns = ["alt", "ref"])

	print(f"{len(CS)} sites loaded.", file = sys.stderr)

	## prefilter poor quality sites

	# 1. excess fraction of MAPQ0 reads at pileup
	frac_mapq0 = CS["mapq0_reads"]/CS["total_reads"] # NOTE: M1 doesnt report sites with cov=0 if not run in forcecalling mode
	mapq_pass_idx = frac_mapq0 <= args.max_frac_mapq0
	print("{} sites with >{}% of MAPQ0 reads will be dropped.".format(len(CS) - mapq_pass_idx.sum(), args.max_frac_mapq0*100), file = sys.stderr)

	# 2. excess fraction of tumor reads pre-filtered by MuTect
	tumor_total_reads = CS["total_reads"] - CS.loc[:, ["n_refcount", "n_altcount"]].sum(1)
	frac_prefiltered = 1 - CS.loc[:, ["t_refcount", "t_altcount"]].sum(1)/tumor_total_reads
	prefilter_pass_idx = frac_prefiltered <= args.max_frac_prefiltered
	print("{} sites with >{}% of prefiltered reads will be dropped.".format(len(CS) - prefilter_pass_idx.sum(), args.max_frac_prefiltered*100), file = sys.stderr)

	# print(CS.loc[~prefilter_pass_idx].head(50), file = sys.stderr)

	# 3. too few reads overall
	tum_cov_idx = (CS["t_altcount"] + CS["t_refcount"] >= args.min_tumor_depth)
	print(f"{(~tum_cov_idx).sum()} sites not sufficiently covered in tumor (cutoff {args.min_tumor_depth}x) will be dropped.", file = sys.stderr)
	print("{} total sites will be dropped; ".format(len(CS) - (mapq_pass_idx & prefilter_pass_idx & tum_cov_idx).sum()), file = sys.stderr, end = "")

	# perform filtering
	CS = CS.loc[mapq_pass_idx & prefilter_pass_idx & tum_cov_idx]

	CS = CS.drop(columns = ["mapq0_reads", "total_reads"])
	print("{} passing sites.".format(CS.shape[0]), file = sys.stderr)

	# load in SNP list
	if args.s is not None:
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
	else:
		H = CS

	good_idx = None
	if not args.use_tonly_genotyper:
		# compute which sites in the SNP list are confidently heterozygous in the normal
		A = H["n_altcount"].values[:, None]
		B = H["n_refcount"].values[:, None]
		# bdens = \int_{af_lb}^{af_ub} df beta(f | n_alt + 1, n_ref + 1)
		H["bdens"] = s.beta.cdf(args.af_ub, A + 1, B + 1) - s.beta.cdf(args.af_lb, A + 1, B + 1)
		# log posterior ratio (alternate method of genotyping; true positive rate is stable WRT coverage)
		H["log_pod"] = np.abs(s.beta.logsf(0.5, A + 1, B + 1) - s.beta.logcdf(0.5, A + 1, B + 1))

		# compute which sites are confidently homozygous alt. in the normal
		H["prob_homalt"] = 1 - s.beta.cdf(0.95, H["n_altcount"].values[:, None] + 1, H["n_refcount"].values[:, None] + 1)

		# save tumor het coverage at good sites to file (GATK GetHetCoverage format)
		if args.use_pod_genotyper:
			good_idx = ( H["log_pod"] < args.log_pod_threshold ) & ( H["n_altcount"]+H["n_refcount"] >= args.pod_min_depth )
		else:
			good_idx = H["bdens"] > args.dens

	# genotype sites based only on the tumor. rather than identifying sites that
	# are confidently homozygous in the normal, we identify sites that are confidently
	# NOT homozygous in the tumor.
	else:
		H["prob_homalt"] = s.beta.sf(0.98, H["t_altcount"].values[:, None] + 1, H["t_refcount"].values[:, None] + 1)
		H["prob_homref"] = s.beta.cdf(0.02, H["t_altcount"].values[:, None] + 1, H["t_refcount"].values[:, None] + 1)

		good_idx = (H["prob_homalt"] < 0.01) & (H["prob_homref"] < 0.1)

	# save tumor het coverage to file
	print("Identified {} high quality het sites in normal.".format(good_idx.sum()), file = sys.stderr)
	H.loc[good_idx, ["chr", "pos", "t_refcount", "t_altcount"]].rename(columns = { "chr" : "CONTIG", "pos" : "POSITION", "t_refcount" : "REF_COUNT", "t_altcount" : "ALT_COUNT" }).to_csv(args.o + ".tumor.tsv", sep = "\t", index = False)

	# save normal het coverage at good sites to file
	H.loc[good_idx, ["chr", "pos", "n_refcount", "n_altcount"]].rename(columns = { "chr" : "CONTIG", "pos" : "POSITION", "n_refcount" : "REF_COUNT", "n_altcount" : "ALT_COUNT" }).to_csv(args.o + ".normal.tsv", sep = "\t", index = False)

	# if requested, save genotype file as TSV (23andme style) 
	if args.g:
		het_idx = good_idx
		hom_idx = (H["prob_homalt"] > args.dens) if not args.use_tonly_genotyper else \
			      (~good_idx & (H["prob_homalt"] > 0.3)) # require minimum coverage of ~17x
		gen_idx = het_idx | hom_idx
		G = H.loc[gen_idx, ["chr", "pos", "allele"]]

		# add genotype info
		alt_ref = np.array(["A", "C", "G", "T"])[np.c_[(G["allele"].values & 0xC) >> 2, G["allele"].values & 3]]
		alt_ref[hom_idx[gen_idx], 1] = alt_ref[hom_idx[gen_idx], 0] 
		G["genotype"] = np.char.add(alt_ref[:, -1], alt_ref[:, 0])

		# restore original contig names
		# XXX: we should probably do this for the coverage files -- how is our
		#      pipeline OK with not doing this?
		G["chr"] = G["chr"].apply(lambda x: contig_list[x - 1])

		# save
		G.drop(columns = ["allele"]).to_csv(args.o + ".genotype.tsv", sep = "\t", index = False)
