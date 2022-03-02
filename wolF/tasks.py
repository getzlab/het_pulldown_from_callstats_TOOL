import wolf
from wolf import Task

def get_het_coverage_from_callstats(
  callstats_file,
  common_snp_list,
  ref_fasta,
  ref_fasta_idx,
  ref_fasta_dict,
  dens_cutoff = None, # if None, automatically infer
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
    if dens_cutoff is None and normal_bam is not None:
        if normal_bai is None:
            raise ValueError("Must explicitly specify BAM index along with BAM!")

        cutoff = Task(
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

    return Task(
        name = "get_het_coverage_from_callstats",
        inputs = {
            "callstats_file" : callstats_file,
            "common_snp_list" : common_snp_list,
            "ref_fasta" : ref_fasta,
            "ref_fasta_idx" : ref_fasta_idx,
            "ref_fasta_dict" : ref_fasta_dict,
            "beta_dens_cutoff" : cutoff
        },
        outputs = {
            "tumor_hets" : "het_coverage.tumor.tsv",
            "normal_hets" : "het_coverage.normal.tsv",
            "normal_genotype" : "het_coverage.genotype.tsv"
        },
        script = "hetpull.py -g -c ${callstats_file} -s ${common_snp_list} -r ${ref_fasta} -o het_coverage --dens ${beta_dens_cutoff}",
        resources = { "mem" : "4G" },
        docker = "gcr.io/broad-getzlab-workflows/het_pulldown_from_callstats:v26"
    )
