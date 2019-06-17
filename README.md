# backmap.pl v0.1

## Description
__Automatic short read mapping and genome size estimate from coverage.__

Automatic mapping of paired and unpaired reads to an assembly with `bwa mem`, execution of `qualimap bamqc` and estimation of genome size from mapped nucleotides divided by mode (>0).
The tools `bwa` and `samtools` need to be in your `$PATH`. The tools `qualimap`, `bedtools` and `Rscript` are optional but needed to create the mapping quality report, coverage histogram and plot of the coverage distribution respectively.

## Dependencies

`backmap.pl` will search for the following executables in your `$PATH`:

Mandatory:
- [bwa (mem)](https://github.com/lh3/bwa): `bwa`
- [samtools](https://github.com/samtools/samtools): `samtools` 

Optional:
- [Qaulimap](http://qualimap.bioinfo.cipf.es/): `qualimap`
- [bedtools](https://bedtools.readthedocs.io/en/latest/) `bedtools`
- [Rscript](https://www.r-project.org/) `Rscript`

## Usage

```
backmap.pl [-a <assembly.fa> {-p <paired_1.fq>,<paired_2.fq> | -u <unpaired.fq>} |
            -b <mapping.bam>]

Mandatory:
	-a STR		Assembly were reads should mapped to in fasta format
	-p STR		Two files with paired reads comma sperated
			Can be specified multiple times
	-u STR		One file with unpaired reads
			Can be specified multiple times
	OR
	-b STR		Bam file to calculate coverage from
			Skips read mapping
			Overrides -nh and -ne

Options: [default]
	-o STR		Output directory [.]
			Will be created if not existing
	-t INT		Number of parallel executed processes [1]
			Affects bwa mem, samtools sort, qualimap bamqc
	-pre STR	Prefix of output files [filename of -a or -b]
	-sort		Sort the bam file (-b) [off]
	-nq		Do not run qualimap bamqc [off]
	-nh		Do not create coverage histogram [off]
			Implies -ne
	-ne		Do not estimate genome size [off]
	-kt		Keep temporary bam files [off]
	-bo STR		Options passed to bwa [-a -c 10000]
			Pass options with -bo "<options>"
	-qo STR		Options passed to qualimap [none]
			Pass options with -qo "<options>"
	-v		Print executed commands to STDERR [off]
	-dry-run	Only print commands to STDERR instead of executing [off]

	-h or -help	Print this help and exit
	-version	Print version number and exit
```

## Citation
If you use this tool please cite the dependencies.
