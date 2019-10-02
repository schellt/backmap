# backmap.pl v0.2

## Description
__Automatic read mapping and genome size estimate from coverage.__

Automatic mapping of paired, unpaired, PacBio and Nanopore reads to an assembly with `bwa mem` or `minimap2`, execution of `qualimap bamqc` and estimation of genome size from mapped nucleotides divided by mode (>0).
The tools `samtools`, `bwa` and/or `minimap2` need to be in your `$PATH`. The tools `qualimap`, `multiqc`, `bedtools` and `Rscript` are optional but needed to create the mapping quality report, coverage histogram and plot of the coverage distribution respectively.

## Dependencies

`backmap.pl` will search for the following executables in your `$PATH`:

Mandatory:
- [Number::FormatEng](https://metacpan.org/pod/Number::FormatEng)
- [samtools](https://github.com/samtools/samtools): `samtools`

Short read mapping:
- [bwa (mem)](https://github.com/lh3/bwa): `bwa`

Long read mapping:
- [minimap2](https://github.com/lh3/minimap2): `minimap2`

Optional:
- [Qaulimap](http://qualimap.bioinfo.cipf.es/): `qualimap`
- [MultiQC](https://multiqc.info/): `multiqc`
- [bedtools](https://bedtools.readthedocs.io/en/latest/) `bedtools`
- [Rscript](https://www.r-project.org/) `Rscript`

## Usage

```
backmap.pl [-a <assembly.fa> {-p <paired_1.fq>,<paired_2.fq> | -u <unpaired.fq>} |
            -pb <pacbio.fq> | -ont <ont.fq> } | -b <mapping.bam>]

Mandatory:
	-a STR		Assembly were reads should mapped to in fasta format
	AND AT LEAST ONE OF
	-p STR		Two files with paired Illumina reads comma sperated
	-u STR		Fastq file with unpaired Illumina reads
	-pb STR		Fasta or fastq file with PacBio reads
	-ont STR	Fasta or fastq file with Nanopore reads
	OR
	-b STR		Bam file to calculate coverage from
			Skips read mapping
			Overrides -nh
			Technologies will recognized correctly if filenames end with
			.pb(.sort).bam or .ont(.sort).bam for PacBio and Nanopore respectively.
			Otherwise they are assumed to be from Illumina.
			
	All mandatory options except of -a can be specified multiple times

Options: [default]
	-o STR		Output directory [.]
			Will be created if not existing
	-t INT		Number of parallel executed processes [1]
			Affects bwa mem, samtools sort, qualimap bamqc
	-pre STR	Prefix of output files if -a is used [filename of -a]
	-sort		Sort the bam file(s) (-b) [off]
	-nq		Do not run qualimap bamqc [off]
	-nh		Do not create coverage histogram [off]
			Implies -ne
	-ne		Do not estimate genome size [off]
	-kt		Keep temporary bam files [off]
	-bo STR		Options passed to bwa [-a -c 10000]
	-mo STR		Options passed to minimap [PacBio: -H -x map-pb; ONT: -x map-ont]
	-qo STR		Options passed to qualimap [none]
	Pass options with quotes e.g. -bo "<options>"
	-v		Print executed commands to STDERR [off]
	-dry-run	Only print commands to STDERR instead of executing [off]

	-h or -help	Print this help and exit
	-version	Print version number and exit
```

## Citation
If you use this tool please cite the dependencies as well.
