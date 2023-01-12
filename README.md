# Tool: het pulldown from callstats

Call normal het sites and get their tumor coverage from a MuTect callstats file.

Normal het sites are called using the posterior odds ratio in a coverage-independent manner. 
The het site allele fraction is assumed to be beta distributed with a flat prior (1, 1):

```
af ~ beta(n_alt + 1, n_ref + 1)
```

A site is considered heterozygous if

```
abs(log(P[af ≥ 0.5] / P[af ≤ 0.5])) < 2.5
```

The threshold can be changed using `--log_pod_threshold`.

Normal het sites can also be called using a simple beta distribution test.
To use this method, specify `--use_beta_density`.
A site is considered heterozygous if

```
Prob[0.4 < af < 0.6] > 0.7
```

The interval bounds (0.4, 0.6) and probability cutoff (0.7) can be overridden.


## Usage:

```
usage: hetpull.py [-h] -c callstats_in -r ref_in -o output_prefix 
                  [-s snplist_in] [-g] 
                  [--use_pod_genotyper] [--log_pod_threshold threshold]
                  [--use_beta_density] [--af_lb lowerbound] [--af_ub upperbound] [--dens cutoff]
                  [--max_frac_mapq0 mapq0_frac] [--max_frac_prefiltered prefiltered_frac]

```

### Required arguments

* `-c callstats_in`: Path to MuTect 1 callstats file
* `-r ref_in`: Path to reference FASTA (directory must also contain FASTA index)
* `-o output_prefix`: Prefix of tumor/normal het site coverage. Tumor coverage will be output to `<prefix>.tumor.tsv`, normal coverage to `<prefix>.normal.tsv`.

### Optional arguments

* `-s snplist_in`: Path to SNP site list, for output of only specific SNP sites. Should be formatted GATK-style, i.e.

   ```
   <chr>    <pos>    <pos>    <strand>    <alt>/<ref>
   ```
   
  Note that the second `<pos>` and `<strand>` can be arbitrary values; all that matters is the number of columns.
  This is to ensure compatibility with GATK's `GetHetCoverage` tool.
* `-g`: Output genotype file (default FALSE)
* `--use_pod_genotyper`: Use posterior odds method for genotyping (default)
* `--log_pod_threshold threshold`: Exclude positions above Log Posterior Odds threshold (default 2.5)
* `--use_beta_density`: Use beta distribution density test for genotyping
* `--af_lb lowerbound`: Lower bound on beta distribution interval (default 0.4)
* `--af_ub upperbound`: Upper bound on beta distribution interval (default 0.6)
* `--dens cutoff`: Probability cutoff (default 0.7)
* `--max_frac_mapq0 mapq0_frac`: Exclude all positions with more than this percentage of MAPQ0 alignments (default 0.05)
* `--max_frac_prefiltered prefiltered_frac`: Exclude all positions with more than this fraction of total reads prefiltered by MuTect (default 0.1)