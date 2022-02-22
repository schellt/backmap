# backmap.pl v0.5

## Description
__Automatic read mapping and genome size estimation from coverage.__

Automatic mapping of paired, unpaired, PacBio and Nanopore reads to an assembly with `bwa mem` or `minimap2`, execution of `qualimap bamqc` and estimation of genome size from mapped nucleotides divided by mode of the coverage distribution (>0). This method was first pulished in Schell et al. (2017). To show high accuracy and reliability of this method throughout the tree of life, Pfenninger et al. (2021) published a study comparing different estimators. Currently, the estimator N<sub>bm</sub>/m (number of back-mapped bases divided by the modal value of the sequencing depth distribution) is implemented in this script only.  
The tools `samtools`, `bwa` and/or `minimap2` need to be in your `$PATH`. The tools `qualimap`, `multiqc`, `bedtools` and `Rscript` are optional but needed to create the mapping quality report, coverage histogram as well as genome size estimation and to plot of the coverage distribution respectively.

## Dependencies

`backmap.pl` needs the following perl modules and will search for executables in your `$PATH`:

Mandatory:
- [Number::FormatEng](https://metacpan.org/pod/Number::FormatEng)
- [Parallel::Loops](https://metacpan.org/pod/Parallel::Loops)
- [samtools](https://github.com/samtools/samtools): `samtools`

Short read mapping:
- [bwa (mem)](https://github.com/lh3/bwa): `bwa`

Long read mapping:
- [minimap2](https://github.com/lh3/minimap2): `minimap2`

Optional:
- [Qualimap](http://qualimap.bioinfo.cipf.es/): `qualimap`
- [MultiQC](https://multiqc.info/): `multiqc`
- [bedtools](https://bedtools.readthedocs.io/en/latest/): `bedtools`
- [Rscript](https://www.r-project.org/): `Rscript`

## Usage

```
backmap.pl [-a <assembly.fa> {-p <paired_1.fq>,<paired_2.fq> | -u <unpaired.fq>} |
            -pb <clr.fq> | -hifi <hifi.fq> | -ont <ont.fq> } | -b <mapping.bam>]

Mandatory:
	-a STR		Assembly were reads should mapped to in fasta format
	AND AT LEAST ONE OF
	-p STR		Two files with paired Illumina reads comma sperated
	-u STR		Fastq file with unpaired Illumina reads
	-pb STR		Fasta or fastq file with PacBio CLR reads
	-hifi STR	Fasta or fastq file with PacBio HiFi reads
	-ont STR	Fasta or fastq file with Nanopore reads
	OR
	-b STR		Bam file to calculate coverage from
			Skips read mapping
			Overrides -nh
			Technologies will recognized correctly if filenames end with
			.pb(.sort).bam, .hifi(.sort).bam or .ont(.sort).bam for PacBio CLR,
			PacBio HiFi and Nanopore respectively. Otherwise they are assumed to
			be from Illumina.
			
	All mandatory options except of -a can be specified multiple times

Options: [default]
	-o STR		Output directory [.]
			Will be created if not existing
	-t INT		Number of parallel executed processes [1]
			Affects bwa mem, samtools sort/index/view/stats, qualimap bamqc
	-pre STR	Prefix of output files if -a is used [filename of -a]
	-sort		Sort the bam file(s) (-b) [off]
	-nq		Do not run qualimap bamqc [off]
	-nh		Do not create coverage histogram [off]
			Implies -ne
	-ne		Do not estimate genome size [off]
	-kt		Keep temporary bam files [off]
	-bo STR		Options passed to bwa [-a -c 10000]
	-mo STR		Options passed to minimap [CLR: -H -x map-pb; HiFi:  minimap<=2.18
			-x asm20 minimap>2.18 -x map-hifi; ONT: -x map-ont]
	-qo STR		Options passed to qualimap [none]
	Pass options with quotes e.g. -bo "<options>"
	-v		Print executed commands to STDERR [off]
	-dry-run	Only print commands to STDERR instead of executing [off]

	-h or -help	Print this help and exit
	-version	Print version number and exit
```

## Citation
Pfenninger M, Schönenbeck P & Schell T (2021). ModEst: Accurate estimation of genome size from next generation sequencing data. _Molecular ecology resources_, 00, 1–11. <https://doi.org/10.1111/1755-0998.13570>  

Schell T, Feldmeyer B, Schmidt H, Greshake B, Tills O et al. (2017). An Annotated Draft Genome for _Radix auricularia_ (Gastropoda, Mollusca). _Genome Biology and Evolution_, 9(3):585–592, <https://doi.org/10.1093/gbe/evx032>

__If you use this tool please cite the dependencies as well:__

- samtools:  
Li H, Handsaker B, Wysoker A, Fennell T, Ruan J et al. (2009). The Sequence Alignment/Map format and SAMtools. _Bioinformatics_, 25(16):2078–2079, <https://doi.org/10.1093/bioinformatics/btp352>
- bwa mem:  
Li H (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. _arXiv preprint arXiv:1303.3997_.
- minimap2:  
Li H (2018). Minimap2: pairwise alignment for nucleotide sequences. _Bioinformatics_, 34:3094–3100, <https://doi.org/10.1093/bioinformatics/bty191>
- Qualimap:  
Okonechnikov K, Conesa A, García-Alcalde F (2016). Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data. _Bioinformatics_, 32(2):292–294, <https://doi.org/10.1093/bioinformatics/btv566>
- MultiQC:  
Ewels P, Magnusson M, Lundin S, Käller M (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. _Bioinformatics_, 32(19):3047–3048, <https://doi.org/10.1093/bioinformatics/btw354>
- bedtools:  
Quinlan AR, Hall IM (2010). BEDTools: a flexible suite of utilities for comparing genomic features. _Bioinformatics_, 26(6):841–842, <https://doi.org/10.1093/bioinformatics/btq033>
- Rscript:  
R Core Team (2021). R: A Language and Environment for Statistical Computing. <http://www.R-project.org/>
