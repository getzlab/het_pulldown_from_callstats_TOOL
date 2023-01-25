import wolf

class get_het_coverage_from_callstats(wolf.Task):
    inputs = {
        "callstats_file" : None,
        "ref_fasta" : None,
        "ref_fasta_idx" : None,
        "ref_fasta_dict" : None,
        "common_snp_list": "",
        "beta_dens_cutoff" : 0.7,
        "log_pod_threshold" : "2.5",
        "max_frac_mapq0" : "0.05",
        "max_frac_prefiltered" : "0.1",
        "use_pod_genotyper" : True,
        "pod_min_depth" : 10
    }
    def script(self):
        return "hetpull.py -g -c ${callstats_file} -r ${ref_fasta} -o het_coverage " \
               "--dens ${beta_dens_cutoff} --max_frac_mapq0 ${max_frac_mapq0} " \
               "--log_pod_threshold ${log_pod_threshold} --pod_min_depth ${pod_min_depth}" + \
               (" --use_pod_genotyper" if self.conf["inputs"]["use_pod_genotyper"] else " --use_beta_density") + \
               (" -s ${common_snp_list}" if self.conf["inputs"]["common_snp_list"] else "")
    outputs = {
        "tumor_hets" : "het_coverage.tumor.tsv",
        "normal_hets" : "het_coverage.normal.tsv",
        "normal_genotype" : "het_coverage.genotype.tsv"
    }
    resources = { "mem" : "4G" }
    docker = "gcr.io/broad-getzlab-workflows/het_pulldown_from_callstats:v42"

class gather_het_coverage(wolf.Task):
    inputs = {
      "tumor_hets",
      "normal_hets",
      "normal_genotype"
    }
    script = """
    cat <(cat $(head -n1 ${normal_genotype}) | head -n1) \
      <(for f in $(cat ${normal_genotype}); do sed 1d $f; done | sort -k1,1V -k2,2n) > normal_genotype.txt
    cat <(cat $(head -n1 ${normal_hets}) | head -n1) \
      <(for f in $(cat ${normal_hets}); do sed 1d $f; done | sort -k1,1V -k2,2n) > normal_hets.txt
    cat <(cat $(head -n1 ${tumor_hets}) | head -n1) \
      <(for f in $(cat ${tumor_hets}); do sed 1d $f; done | sort -k1,1V -k2,2n) > tumor_hets.txt
    """
    outputs = {
      "tumor_hets" : "tumor_hets.txt",
      "normal_hets" : "normal_hets.txt",
      "normal_genotype" : "normal_genotype.txt",
    }

def get_het_coverage_from_callstats_legacy_workflow(
  callstats_file,
  common_snp_list,
  ref_fasta,
  ref_fasta_idx,
  ref_fasta_dict,
  dens_cutoff = None, # if None, automatically infer
  max_frac_mapq0 = 0.05,
  use_pod_genotyper = False,
  log_pod_threshold = 2.5,
  pod_min_depth = 10,
  normal_bam = None,
  normal_bai = None,
):
    """
    Get het site coverage from a MuTect callstats file
    -- callstats_files: path to MuTect callstats file
    -- common_snp_list: path to list of germline SNPs that will be surveyed
    -- refFasta: path to reference FASTA File
    -- dens_cutoff: beta density cutoff for genotyping normal het sites
    -- normal_bam: path to normal BAM to automatically determine good cutoff based on coverage
    """

    # infer density cutoff from BAM coverage
    if not use_pod_genotyper and dens_cutoff is None and normal_bam is not None:
        if normal_bai is None:
            raise ValueError("Must explicitly specify BAM index along with BAM!")

        cutoff = wolf.Task(
          name = "set_het_density_cutoff",
          inputs = {
            "normal_bam" : normal_bam,
            "normal_bai" : normal_bai
          },
          script = """
# get read length
READ_LEN=$(samtools view ${normal_bam} | head -n100 | cut -f6 | \
perl -ne '@ciglens = split(/[^0-9]+/); $tot = 0;
foreach(@ciglens) { $tot += $_; }
$maxtot = $tot > $maxtot ? $tot : $maxtot;
END { print $maxtot }')
# get approximate mean coverage
COV=$(samtools idxstats ${normal_bam} | awk -v read_len=$READ_LEN '
NR <= 24 { genome_len += $2; n_reads += $3 }
END { print read_len*n_reads/genome_len }
')
# get appropriate density cutoff
python - $COV <<EOF > cutoff.txt
import numpy as np
import scipy.stats as s
import sys
cov = float(sys.argv[1])
print(0.9*(1 - np.diff(s.beta.cdf(np.r_[0.05, 0.4, 0.6, 0.95], cov/2 + 1, cov/2 + 1))[[0, -1]].sum()))
EOF
""",
          outputs = {
            "cutoff" : ("cutoff.txt", wolf.read_file)
          },
          docker = "gcr.io/broad-getzlab-workflows/base_image:v0.0.5"
        )["cutoff"]

    # set to default value of 0.7
    elif dens_cutoff is None:
        cutoff = 0.7

    # use user-specified value
    else:
        if dens_cutoff > 1 or dens_cutoff < 0:
            raise ValueError("Density cutoff must be between 0 and 1!")
        cutoff = dens_cutoff

    return get_het_coverage_from_callstats(
        inputs = {
            "callstats_file" : callstats_file,
            "common_snp_list" : common_snp_list,
            "ref_fasta" : ref_fasta,
            "ref_fasta_idx" : ref_fasta_idx,
            "ref_fasta_dict" : ref_fasta_dict,
            "beta_dens_cutoff" : cutoff,
            "log_pod_threshold" : log_pod_threshold,
            "pod_min_depth" : pod_min_depth,
            "max_frac_mapq0" : max_frac_mapq0
        }
    )
