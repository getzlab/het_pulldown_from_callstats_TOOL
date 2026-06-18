#!/usr/bin/env python

import argparse
import numpy as np
import os
import pandas as pd
import scipy.stats as s
import subprocess
import sys
from capy import seq

from hetmodels import run_snp_mixture_model

REF = 'ref_allele'
ALT = 'alt_allele'
CHROMOSOME = 'contig'
POSITION = 'position'
TOTAL_READS = 'total_reads'
MAP_Q0_READS = 'map_Q0_reads'
T_REF_COUNT = 't_ref_count'
T_ALT_COUNT = 't_alt_count'
N_REF_COUNT = 'n_ref_count'
N_ALT_COUNT = 'n_alt_count'
ALLELE = 'allele'

def parse_args():
	# parse args
	parser = argparse.ArgumentParser(description = "Get het site coverage from MuTect 1 callstats file")
	parser.add_argument("-c", "--callstats", required = True, help = "Path to callstats file", metavar = "callstats_in")
	parser.add_argument("-r", "--ref-fasta", required = True, help = "Path to reference FASTA (directory must contain FASTA index)", metavar = "ref_in")
	parser.add_argument("-o", "--out-prefix", required = True, help = "Het coverage file prefix ('tumor'/'normal' appended)", metavar = "output_prefix")
	parser.add_argument("-s", "--snp-list", help="Path to GATK-formatted SNP site file", metavar="snplist_in")
	parser.add_argument("-g", "--genotype", help = "Output genotype file", action = "store_true")
	parser.add_argument("--mutect", help = "Marks that the input was produced by MuTect", action = "store_true")

	genotyper_parser = parser.add_mutually_exclusive_group(required=False)
	genotyper_parser.add_argument("-m", dest="method", help = "Selection method to use: mixture_model, pod, beta_density", choices=['mixture_model', 'pod', 'beta_density'], metavar="method")
	genotyper_parser.add_argument("--use_pod_genotyper", dest="use_pod_genotyper",help = "(Deprecated use --method argument) Use posterior odds method for genotyping", action = "store_true")
	genotyper_parser.add_argument("--use_beta_density", dest="use_pod_genotyper", help = "(Deprecated use --method argument) Use beta distribution density for genotyping", action = "store_false")
	
	parser.add_argument("--use_tonly_genotyper", help = "Genotype in single sample mode", action = "store_true")
	parser.add_argument("--pod_min_depth", type=int, default=10,
						help="(Deprecated, use min_normal_depth) Any position with total normal coverage below this threshold will not be considered for genotyping")
	parser.add_argument("--min_normal_depth", type=int, default=10,
						help="Any position with total normal coverage below this threshold will not be considered for genotyping")
	parser.add_argument("--min_tumor_depth", type=int, default=1,
						help="Any position with total tumor coverage below this threshold will be discarded")

	parser.add_argument("--log_pod_threshold", type = float, metavar='threshold', default = 2.5)
	parser.add_argument("--af_lb", help = "Lower bound on beta distribution AF interval", default = 0.4, type = float, metavar = "lowerbound")
	parser.add_argument("--af_ub", help = "Upper bound on beta distribution AF interval", default = 0.6, type = float, metavar = "upperbound")
	parser.add_argument("--dens", help = "Beta distribution density threshold to consider a site heterozygous in the normal.", default = 0.7, type = float, metavar = "cutoff")

	parser.add_argument("--max_frac_mapq0", help = "Any position from callstats with more than this percentage of MAPQ0 reads will be excluded from het coverage analysis.", default = 0.05, type = float, metavar = "mapq0_frac")
	parser.add_argument("--max_frac_prefiltered", help = "Any site with more than this fraction of total reads prefiltered by MuTect will be excluded.", default = 0.10, type = float, metavar = "prefiltered_frac")

	parser.set_defaults(use_pod_genotyper=False)
	parser.set_defaults(use_beta_density=False)
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

	if args.use_pod_genotyper or args.use_beta_density:
		raise ValueError("--use_pod_genotyper and --use_beta_density deprecated, use -m argument")

	if not os.path.exists(args.callstats):
		raise FileNotFoundError("Callstats file not found!")
	if args.snp_list is None:
		print("WARNING: without an input SNP site file, all valid het sites will be returned.", file = sys.stderr)
	elif not os.path.exists(args.snp_list):
		raise FileNotFoundError("SNP site file not found!")
	if not os.path.exists(args.ref_fasta):
		raise FileNotFoundError("Reference fasta file not found!")
	if not os.path.exists(args.ref_fasta + '.fai'):
		raise FileNotFoundError("Reference fasta index file not found! (Must be <reference.fa>.fai)")

	return args

def hash_altref(DF):
	return (DF.replace(dict(zip(list("ACGT"), range(0, 4))))@[4, 1]).astype(np.uint8)


def load_callstats_file(cs_file: str, ref_file: str, is_mutect: bool):
	print("Loading callstats file ...", file=sys.stderr)
	if is_mutect:
		# trim callstats (faster to do this on the shell)
		callstats_trimmed = subprocess.Popen("sed '1,2d' {} | cut -f1,2,4,5,16,17,26,27,38,39".format(cs_file), shell=True, stdout=subprocess.PIPE)
		assert callstats_trimmed.stdout is not None, 'Failed to trim MuTect output!'
		# load in callstats
		CS = pd.read_csv(callstats_trimmed.stdout, sep="\t",
						names=[CHROMOSOME, POSITION, REF, ALT, TOTAL_READS, MAP_Q0_READS, T_REF_COUNT, T_ALT_COUNT,
								N_REF_COUNT, N_ALT_COUNT],
						dtype={CHROMOSOME: str, POSITION: np.uint32, TOTAL_READS: np.uint32, MAP_Q0_READS: np.uint32,
								T_REF_COUNT: np.uint32, T_ALT_COUNT: np.uint32, N_REF_COUNT: np.uint32,
								N_ALT_COUNT: np.uint32}
						)
	else:
		CS = pd.read_csv(cs_file, sep='\t', comment='#')

	contig_list = pd.read_csv(ref_file + '.fai', sep='\t', usecols=[0], names=["contig"])["contig"].tolist()
	CS[CHROMOSOME] = CS[CHROMOSOME].apply(lambda x: contig_list.index(x) + 1).astype(np.uint8)
	CS["gpos"] = seq.chrpos2gpos(CS[CHROMOSOME], CS[POSITION], ref=ref_file)
	CS[ALLELE] = hash_altref(CS.loc[:, [ALT, REF]])
	CS = CS.drop(columns=[ALT, REF])

	return CS


def apply_prefilters(CS,max_frac_mapq0,max_frac_prefiltered,min_tumor_depth):
	mask = np.full(len(CS), True)
	
	# 1. excess fraction of MAPQ0 reads at pileup
	if TOTAL_READS in CS and MAP_Q0_READS in CS:
		frac_mapq0 = CS[MAP_Q0_READS]/CS[TOTAL_READS] # NOTE: M1 doesnt report sites with cov=0 if not run in forcecalling mode
		mapq_pass_idx = frac_mapq0 <= max_frac_mapq0
		mask &= mapq_pass_idx
		print("{} sites with >{}% of MAPQ0 reads will be dropped.".format(len(CS) - mapq_pass_idx.sum(), max_frac_mapq0*100), file = sys.stderr)

	# 2. excess fraction of tumor reads pre-filtered by MuTect
	if TOTAL_READS in CS and N_REF_COUNT in CS and N_ALT_COUNT in CS and T_REF_COUNT in CS and T_ALT_COUNT in CS:
		tumor_total_reads = CS[TOTAL_READS] - CS.loc[:, [N_REF_COUNT, N_ALT_COUNT]].sum(1)
		frac_prefiltered = 1 - CS.loc[:, [T_REF_COUNT, T_ALT_COUNT]].sum(1)/tumor_total_reads
		prefilter_pass_idx = frac_prefiltered <= max_frac_prefiltered
		mask &= prefilter_pass_idx
		print("{} sites with >{}% of prefiltered reads will be dropped.".format(len(CS) - prefilter_pass_idx.sum(), max_frac_prefiltered*100), file = sys.stderr)


	# 3. too few reads overall
	if T_REF_COUNT in CS and T_ALT_COUNT in CS:
		tum_cov_idx = (CS[T_ALT_COUNT] + CS[T_REF_COUNT] >= min_tumor_depth)
		mask &= tum_cov_idx
		print(f"{(~tum_cov_idx).sum()} sites not sufficiently covered in tumor (cutoff {min_tumor_depth}x) will be dropped.", file = sys.stderr)

	print("{} total sites will be dropped; ".format(len(CS) - mask.sum()), file = sys.stderr, end = "")
	# perform filtering
	CS = CS.loc[mask]

	CS = CS.drop(columns = [MAP_Q0_READS, TOTAL_READS])

	return(CS)

if __name__ == "__main__":
	args = parse_args()

	contig_list = pd.read_csv(args.ref_fasta + '.fai', sep='\t', usecols = [0], names=["contig"])["contig"].tolist()
	CS = load_callstats_file(args.callstats, args.ref_fasta, args.mutect)

	print(f"{len(CS)} sites loaded.", file = sys.stderr)

	## prefilter poor quality sites

	CS = apply_prefilters(CS,max_frac_mapq0=args.max_frac_mapq0,
						  max_frac_prefiltered=args.max_frac_prefiltered,
						  min_tumor_depth=args.min_tumor_depth)

	print("{} passing sites.".format(CS.shape[0]), file = sys.stderr)

	# load in SNP list
	if args.snp_list is not None:
		print("Loading SNP list ...", file = sys.stderr)
		H = pd.read_csv(args.snp_list, sep = "\t", comment = "@",
		  names = [CHROMOSOME, POSITION, "x", "y", ALLELE],
		  dtype = { CHROMOSOME : str, POSITION : np.uint32, "x" : np.uint32, "y" : str, ALLELE : str },
		).drop(columns = ["x", "y"])
		H[CHROMOSOME] = H[CHROMOSOME].apply(lambda x: contig_list.index(x) + 1).astype(np.uint8)
		H["gpos"] = seq.chrpos2gpos(H[CHROMOSOME], H[POSITION], ref = args.ref_fasta)
		H[ALLELE] = hash_altref(H[ALLELE].str.extract(r"(.)/(.)"))
		print("{} sites loaded.".format(H.shape[0]), file = sys.stderr)

		# merge
		print("Pulling down SNP site coverage from callstats ...", file = sys.stderr)
		H = H.merge(CS.drop(columns = [CHROMOSOME, POSITION]), left_on = ["gpos", ALLELE], right_on = ["gpos", ALLELE], how = "inner")
		print("{} covered SNP sites identified.".format(H.shape[0]), file = sys.stderr)
	else:
		H = CS

	good_idx = None

	if not args.use_tonly_genotyper:
		A = H[N_ALT_COUNT].values[:, None]
		B = H[N_REF_COUNT].values[:, None]

	else:
		A = H[T_ALT_COUNT].values[:, None]
		B = H[T_REF_COUNT].values[:, None]
		print('Using Tumor Only Genotyping')

	if args.method == "mixture_model":
		outs = run_snp_mixture_model(B,A)


		if args.use_tonly_genotyper:
			outs = run_snp_mixture_model(B,A)
			H[outs['snp_prob'].columns] = outs['snp_prob'].values
			good_idx = (H['prob_het'] + H['prob_other']) > .9999
		else:
			outs = run_snp_mixture_model(B,A,include_noise_component=True,fix_noise_component=1e-5)
			H[outs['snp_prob'].columns] = outs['snp_prob'].values
			good_idx = H['prob_het'] > .999
            
	elif not args.use_tonly_genotyper:
		# compute which sites in the SNP list are confidently heterozygous in the normal
		A = H[N_ALT_COUNT].values[:, None]
		B = H[N_REF_COUNT].values[:, None]
		# bdens = \int_{af_lb}^{af_ub} df beta(f | n_alt + 1, n_ref + 1)
		H["bdens"] = s.beta.cdf(args.af_ub, A + 1, B + 1) - s.beta.cdf(args.af_lb, A + 1, B + 1)
		# log posterior ratio (alternate method of genotyping; true positive rate is stable WRT coverage)
		H["log_pod"] = np.abs(s.beta.logsf(0.5, A + 1, B + 1) - s.beta.logcdf(0.5, A + 1, B + 1))

		# compute which sites are confidently homozygous alt. in the normal
		H["prob_homalt"] = 1 - s.beta.cdf(0.95, H[N_ALT_COUNT].values[:, None] + 1, H[N_REF_COUNT].values[:, None] + 1)

		# save tumor het coverage at good sites to file (GATK GetHetCoverage format)
		if args.method=='pod':
			good_idx = ( H["log_pod"] < args.log_pod_threshold ) & ( H[N_ALT_COUNT]+H[N_REF_COUNT] >= args.pod_min_depth )
		else:
			good_idx = H["bdens"] > args.dens

	# genotype sites based only on the tumor. rather than identifying sites that
	# are confidently homozygous in the normal, we identify sites that are confidently
	# NOT homozygous in the tumor.
	else:
		H["prob_homalt"] = s.beta.sf(0.98, H[T_ALT_COUNT].values[:, None] + 1, H[T_REF_COUNT].values[:, None] + 1)
		H["prob_homref"] = s.beta.cdf(0.02, H[T_ALT_COUNT].values[:, None] + 1, H[T_REF_COUNT].values[:, None] + 1)

		good_idx = (H["prob_homalt"] < 0.01) & (H["prob_homref"] < 0.1)

	# save all possible het sites to file
	H.to_csv(args.out_prefix + ".all_sites.tsv", sep = "\t", index = False)

	# save tumor het coverage to file
	print("Identified {} high quality het sites in normal.".format(good_idx.sum()), file = sys.stderr)
	H.loc[good_idx, [CHROMOSOME, POSITION, T_REF_COUNT, T_ALT_COUNT]].rename(columns = { CHROMOSOME : "CONTIG", POSITION : "POSITION", T_REF_COUNT : "REF_COUNT", T_ALT_COUNT : "ALT_COUNT" }).to_csv(args.out_prefix + ".tumor.tsv", sep = "\t", index = False)

	# save normal het coverage at good sites to file
	H.loc[good_idx, [CHROMOSOME, POSITION, N_REF_COUNT, N_ALT_COUNT]].rename(columns = { CHROMOSOME : "CONTIG", POSITION : "POSITION", N_REF_COUNT : "REF_COUNT", N_ALT_COUNT : "ALT_COUNT" }).to_csv(args.out_prefix + ".normal.tsv", sep = "\t", index = False)

	# if requested, save genotype file as TSV (23andme style) 
	if args.genotype:
		het_idx = good_idx
		hom_idx = (H["prob_homalt"] > args.dens) if not args.use_tonly_genotyper else \
			      (~good_idx & (H["prob_homalt"] > 0.3)) # require minimum coverage of ~17x
		gen_idx = het_idx | hom_idx
		G = H.loc[gen_idx, [CHROMOSOME, POSITION, ALLELE]]

		# add genotype info
		alt_ref = np.array(["A", "C", "G", "T"])[np.c_[(G[ALLELE].values & 0xC) >> 2, G[ALLELE].values & 3]]
		alt_ref[hom_idx[gen_idx], 1] = alt_ref[hom_idx[gen_idx], 0] 
		G["genotype"] = np.char.add(alt_ref[:, -1], alt_ref[:, 0])

		# restore original contig names
		# XXX: we should probably do this for the coverage files -- how is our
		#      pipeline OK with not doing this?
		G[CHROMOSOME] = G[CHROMOSOME].apply(lambda x: contig_list[x - 1])

		# save
		G.drop(columns = [ALLELE]).to_csv(args.out_prefix + ".genotype.tsv", sep = "\t", index = False)

