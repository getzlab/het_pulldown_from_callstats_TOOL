# Method: het pulldown from callstats

Call normal het sites and get their tumor coverage from a MuTect callstats file.

Normal het sites are called using a simple beta distribution test. 
The het site allele fraction is assumed to be beta distributed with a flat prior (1, 1):

```
af ~ beta(n_alt + 1, n_ref + 1)
```

A site is considered heterozygous if

```
Prob[0.4 < af < 0.6] > 0.7
```

The interval bounds (0.4, 0.6) and probability cutoff (0.7) can be overridden.

## Usage:

```
hetpull.py [-h] -c callstats_in -s snplist_in [-o output_path]
           [--af_lb lowerbound] [--af_ub upperbound] [--dens cutoff]
```

### Required arguments

* `-c callstats_in`: Path to MuTect 1 callstats file
* `-s snplist_in`: Path to SNP site list. Should be formatted GATK-style, i.e.

   ```
   <chr>    <pos>    <pos>    <strand>    <alt>/<ref>
   ```
   
Note that the second `<pos>` and `<strand>` can be arbitrary values; all that matters is the number of columns.
This is to ensure compatibility with GATK's `GetHetCoverage` tool.

### Optional arguments

* `-o output_path`: Path to coverage output. Will default to `stdout` if not specified.
* `--af_lb lowerbound`: Lower bound on beta distribution interval (default 0.4)
* `--af_ub upperbound`: Upper bound on beta distribution interval (default 0.6)
* `--dens cutoff`: Probability cutoff (default 0.7)
