#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';
use IPC::Cmd qw[can_run run];
use Number::FormatEng qw(:all);
use Parallel::Loops;

my $version = "0.3";

sub print_help{
	print STDOUT "\n";
	print STDOUT "backmap.pl v$version\n";
	print STDOUT "\n";
	print STDOUT "Description:\n";
	print STDOUT "\tAutomatic mapping of paired, unpaired, PacBio and Nanopore reads to an\n\tassembly, execution of qualimap bamqc, multiqc and estimation of genome size\n\tfrom mapped nucleotides and peak coverage.\n\tThe tools bwa, minimap2, samtools, qualimap, multiqc, bedtools and Rscript need to be\n\tin your \$PATH.\n";
	print STDOUT "\n";
	print STDOUT "Usage:\n";
	print STDOUT "\tbackmap.pl [-a <assembly.fa> {-p <paired_1.fq>,<paired_2.fq> | -u <unpaired.fq> |\n";
	print STDOUT "\t            -pb <pacbio.fq> | -ont <ont.fq> } | -b <mapping.bam>]\n";
	print STDOUT "\n";
	print STDOUT "Mandatory:\n";
	print STDOUT "\t-a STR\t\tAssembly were reads should mapped to in fasta format\n";
	print STDOUT "\tAND AT LEAST ONE OF\n";
	print STDOUT "\t-p STR\t\tTwo fastq files with paired Illumina reads comma sperated\n";
	print STDOUT "\t-u STR\t\tFastq file with unpaired Illumina reads\n";
	print STDOUT "\t-pb STR\t\tFasta or fastq file with PacBio reads\n";
	print STDOUT "\t-ont STR\tFasta or fastq file with Nanopore reads\n";
	print STDOUT "\tOR\n";
	print STDOUT "\t-b STR\t\tBam file to calculate coverage from\n";
	print STDOUT "\t\t\tSkips read mapping\n";
	print STDOUT "\t\t\tOverrides -nh\n";
	print STDOUT "\t\t\tTechnologies will recognized correctly if filenames end with\n\t\t\t.pb(.sort).bam or .ont(.sort).bam for PacBio and Nanopore respectively.\n\t\t\tOtherwise they are assumed to be from Illumina.\n";
	print STDOUT "\tAll mandatory options except of -a can be specified multiple times\n";
	print STDOUT "\n";
	print STDOUT "Options: [default]\n";
	print STDOUT "\t-o STR\t\tOutput directory [.]\n";
	print STDOUT "\t\t\tWill be created if not existing\n";
	print STDOUT "\t-t INT\t\tNumber of parallel executed processes [1]\n";
	print STDOUT "\t\t\tAffects bwa mem, samtools sort, qualimap bamqc\n";
	print STDOUT "\t-pre STR\tPrefix of output files if -a is used [filename of -a]\n";
	print STDOUT "\t-sort\t\tSort the bam file(s) (-b) [off]\n";
	print STDOUT "\t-nq\t\tDo not run qualimap bamqc [off]\n";
	print STDOUT "\t-nh\t\tDo not create coverage histogram [off]\n";
	print STDOUT "\t\t\tImplies -ne\n";
	print STDOUT "\t-ne\t\tDo not estimate genome size [off]\n";
	print STDOUT "\t-kt\t\tKeep temporary bam files [off]\n";
	print STDOUT "\t-bo STR\t\tOptions passed to bwa [-a -c 10000]\n";
	print STDOUT "\t-mo STR\t\tOptions passed to minimap [PacBio: -H -x map-pb; ONT: -x map-ont]\n";
	print STDOUT "\t-qo STR\t\tOptions passed to qualimap [none]\n";
	print STDOUT "\tPass options with quotes e.g. -bo \"<options>\"\n";
	print STDOUT "\t-v\t\tPrint executed commands to STDERR [off]\n";
	print STDOUT "\t-dry-run\tOnly print commands to STDERR instead of executing [off]\n";
	print STDOUT "\n";
	print STDOUT "\t-h or -help\tPrint this help and exit\n";
	print STDOUT "\t-version\tPrint version number and exit\n";
	exit;
}

sub exe_cmd{
	my ($cmd,$verbose,$dry) = @_;
	if($verbose == 1){
		print STDERR "CMD\t$cmd\n";
	}
	if($dry == 0){
		system("$cmd") == 0 or die "ERROR\tsystem $cmd failed: $?";
	}
}

sub round_format_pref{
	my ($number) = @_;
	$number = format_pref($number);
	my $format_number = substr($number,0,-1);
	$format_number = sprintf("%.2f", $format_number);
	my $pref = substr($number,-1);
	my $return = $format_number . $pref;
	return $return;
}

my $out_dir = abs_path("./");
my $assembly_path = "";
my $assembly = "";
my @paired = ();
my @unpaired = ();
my @pb = ();
my @ont = ();
my $threads = 1;
my $prefix = "";
my $verbose = 0;
my $bwa_opts = "-a -c 10000 ";
my $minimap_opts = "";
my $qm_opts = "";
my $create_histo_switch = 1;
my $estimate_genome_size_switch = 1;
my $run_bamqc_switch = 1;
my $keep_tmp = 0;
my $dry = 0;
my @bam = ();
my $sort_bam_switch = 0;
my $cmd;
my $mkdir_cmd;
my $home = `echo \$HOME`;
chomp $home;

my $input_error = 0;

if(scalar(@ARGV==0)){
	print_help;
}

for (my $i = 0; $i < scalar(@ARGV);$i++){
	if ($ARGV[$i] eq "-o"){
		$out_dir = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-a"){
		if($assembly ne ""){
			print STDERR "ERROR\tSpecify -a just once\n";
			$input_error = 1;
		}
		$assembly = (split /\//,$ARGV[$i+1])[-1];
		$assembly_path = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-p"){
		push(@paired,$ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-u"){
		push(@unpaired,$ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-pb"){
		push(@pb,$ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-ont"){
		push(@ont,$ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-b"){
		push(@bam,$ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-sort"){
		$sort_bam_switch = 1;
	}
	if ($ARGV[$i] eq "-t"){
		$threads = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-pre"){
		$prefix = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-v"){
		$verbose = 1;
	}
	if ($ARGV[$i] eq "-bo"){
		$bwa_opts = $ARGV[$i+1] . " ";	#nonsense flags are skipped from bwa
		$ARGV[$i+1] = "\'$ARGV[$i+1]\'";
	}
	if ($ARGV[$i] eq "-mo"){
		$minimap_opts = $ARGV[$i+1] . " ";
		$ARGV[$i+1] = "\'$ARGV[$i+1]\'";
	}
	if ($ARGV[$i] eq "-qo"){
		$qm_opts = $ARGV[$i+1] . " ";	#nonsense flags are skipped from qualimap
		$ARGV[$i+1] = "\'$ARGV[$i+1]\'";
	}
	if ($ARGV[$i] eq "-nq"){
		$run_bamqc_switch = 0;
	}
	if ($ARGV[$i] eq "-nh"){
		$create_histo_switch = 0;
		$estimate_genome_size_switch = 0;
	}
	if ($ARGV[$i] eq "-ne"){
		$estimate_genome_size_switch = 0;
	}
	if($ARGV[$i] eq "-kt"){
		$keep_tmp = 1;
	}
	if ($ARGV[$i] eq "-dry-run"){
		$dry = 1;
		$verbose = 1;
	}
	if ($ARGV[$i] eq "-h" or $ARGV[$i] eq "-help"){
		print_help;
	}
	if ($ARGV[$i] eq "-version"){
		print STDERR $version . "\n";
		exit;
	}
}

print STDERR "CMD\t" . $0 . " " . join(" ",@ARGV) . "\n";

if($assembly_path ne "" and scalar(@bam) > 0){
	print STDERR "ERROR\tSpecify either -a or -b\n";
	$input_error = 1;
}

if($assembly_path eq "" and scalar(@bam) == 0){
	print STDERR "ERROR\tSpecify either -a or -b\n";
	$input_error = 1;
}

if($assembly_path ne "" and scalar(@bam) == 0){
	if(scalar(@paired) > 0 or scalar(@unpaired) > 0){
		if(not defined(can_run("bwa"))){
			print STDERR "ERROR\tbwa is not in your \$PATH\n";
			$input_error = 1;
		}
	}
	if(scalar(@pb) > 0 or scalar(@ont) > 0){
		if(not defined(can_run("minimap2"))){
			print STDERR "ERROR\tminimap2 is not in your \$PATH\n";
			$input_error = 1;
		}
	}
	if(not defined(can_run("samtools"))){
		print STDERR "ERROR\tsamtools is not in your \$PATH\n";
		$input_error = 1;
	}
}

if(not defined(can_run("qualimap")) and $run_bamqc_switch == 1){
	print STDERR "INFO\tqualimap is not in your \$PATH and will not be executed\n";
}

if(not defined(can_run("bedtools"))){
	print STDERR "INFO\tbedtools is not in your \$PATH\n";
	if($estimate_genome_size_switch == 1){
		print STDERR "WARNING\tGenome size estimation not possible\n";
		$estimate_genome_size_switch = 0;
	}
	if(scalar(@bam) > 0){
		print STDERR "ERROR\tGenome size estimation not possible\n";
		$input_error = 1;
	}
}

if(not defined(can_run("Rscript"))){
	print STDERR "INFO\tRscript is not in your \$PATH\n";
	if($create_histo_switch == 1){
		print STDERR "WARNING\tPlotting not possible\n";
	}
}

if(-f "$out_dir"){
	print STDERR "ERROR\tOutput directory $out_dir is already a file!\n";
	$input_error = 1;
}

if($assembly_path ne ""){
	if(not -f $assembly_path){
		print STDERR "ERROR\tFile $assembly_path does not exist!\n";
		$input_error = 1;
	}
	if(scalar(@paired) == 0 and scalar(@unpaired) == 0 and scalar(@pb) == 0 and scalar(@ont) == 0){
		print STDERR "ERROR\tNo reads specified!\n";
		$input_error = 1;
	}
}

if($threads !~ m/^\d+$/ or $threads < 1){
	print STDERR "ERROR\tThreads is no integer >= 1!\n";
	$input_error = 1;
}

if ($input_error == 1){
	print STDERR "ERROR\tInput error detected!\n";
	exit 1;
}

my $samtools_threads = $threads - 1;

if(not -d "$out_dir"){
	print STDERR "INFO\tCreating output directory $out_dir\n";
	$mkdir_cmd="mkdir -p $out_dir";
}

if($prefix eq ""){
	if($assembly ne ""){
		$prefix = $assembly;
		print STDERR "INFO\tSetting prefix to $prefix\n";
	}
	else{
		print STDERR "INFO\tSetting prefix to name(s) of corresponding bam\n";
	}
}

if(scalar(@bam) > 0){
	if(defined(can_run("bedtools"))){
		if($create_histo_switch == 0 or $estimate_genome_size_switch == 0){
			$create_histo_switch = 1;
			print STDERR "INFO\tCreating coverage histogram\n";
		}
	}
}

my %paired_filter;

if($assembly_path ne ""){
	foreach(@paired){
		my @pair = split(/,/,$_);
		foreach(@pair){
			if($_ =~ m/^~/){
				$_ =~ s/^~/$home/;      #~ is translated by bash into $HOME. This does not work if there is no space infront. That means if the second file starts with "~" it will not be recognized even though it exists
			}
		}
		if(scalar(@pair) != 2){
			print STDERR "INFO\tNot a pair: $_ - skipping these file(s)\n";
		}
		else{
			my $file_error = 0;
			if(not -f "$pair[0]"){
				print STDERR "INFO\tNo file $pair[0] - skipping pair $_\n";
				$file_error = 1;
			}
			if(not -f "$pair[1]"){
				print STDERR "INFO\tNo file $pair[1] - skipping pair $_\n";
				$file_error = 1;
			}
			if($file_error == 0){
				if(exists($paired_filter{abs_path($pair[0]) . "," . abs_path($pair[1])})){
					print STDERR "INFO\tPair " . abs_path($pair[0]) . "," . abs_path($pair[1]) . " already specified\n";
				}
				else{
					$paired_filter{abs_path($pair[0]) . "," . abs_path($pair[1])} = 1;
				}
			}
		}
	}
}

my %unpaired_filter;

if($assembly_path ne ""){
	foreach(@unpaired){
		if(not -f "$_"){
			print STDERR "INFO\tNo file $_ - skipping this file\n";
		}
		else{
			if(exists($unpaired_filter{abs_path($_)})){
				print STDERR "INFO\tFile " . abs_path($_) . " already specified\n";
			}
			else{
				$unpaired_filter{abs_path($_)} = 1;
			}
		}
	}
}

my %pb_filter;

if($assembly_path ne ""){
	foreach(@pb){
		if(not -f "$_"){
			print STDERR "INFO\tNo file $_ - skipping this file\n";
		}
		else{
			if(exists($pb_filter{abs_path($_)})){
				print STDERR "INFO\tFile " . abs_path($_) . " already specified\n";
			}
			else{
				$pb_filter{abs_path($_)} = 1;
			}
		}
	}
}

my %ont_filter;

if($assembly_path ne ""){
	foreach(@ont){
		if(not -f "$_"){
			print STDERR "INFO\tNo file $_ - skipping this file\n";
		}
		else{
			if(exists($ont_filter{abs_path($_)})){
				print STDERR "INFO\tFile " . abs_path($_) . " already specified\n";
			}
			else{
				$ont_filter{abs_path($_)} = 1;
			}
		}
	}
}

my %bam_filter;

if(scalar(@bam) > 0){
	foreach(@bam){
		if(not -f "$_"){
			print STDERR "INFO\tNo file $_ - skipping this file\n";
		}
		else{
			if(exists($bam_filter{abs_path($_)})){
				print STDERR "INFO\tFile " . abs_path($_) . " already specified\n";
			}
			else{
				$bam_filter{abs_path($_)} = 1;
			}
		}
	}
}

if($assembly_path ne ""){
	if(scalar(keys(%paired_filter)) == 0 and scalar(keys(%unpaired_filter)) == 0 and scalar(keys(%pb_filter)) == 0 and scalar(keys(%ont_filter)) == 0){
		print STDERR "ERROR\tNo existing read files specified!\n";
		exit 1;
	}
}
else{
	if(scalar(keys(%bam_filter)) == 0){
		print STDERR "ERROR\tNo existing bam files specified!\n";
		exit 1;
	}
}

my $index_path = $out_dir . "/" . $assembly;

if($assembly_path ne ""){
	if(scalar(keys(%paired_filter)) > 0 or scalar(keys(%unpaired_filter)) > 0){
		if(not -f $assembly_path . ".amb" or not -f $assembly_path . ".ann" or not -f $assembly_path . ".bwt" or not -f $assembly_path . ".pac" or not -f $assembly_path . ".sa"){
			$cmd = "bwa index -p $index_path $assembly_path > $out_dir/$prefix\_bwa_index.log 2> $out_dir/$prefix\_bwa_index.err";
		}
		else{
			print STDERR "INFO\tIndex files already existing\n";
			$index_path = $assembly_path;
		}
	}
}

my $bwa_version;
if(not defined(can_run("bwa"))){
	$bwa_version = "not detected";
}
else{
	$bwa_version = `bwa 2>&1 | head -3 | tail -1 | sed 's/^Version: //'`;
	chomp $bwa_version;
}

my $minimap_version;
if(not defined(can_run("minimap2"))){
	$minimap_version = "not detected";
}
else{
	$minimap_version = `minimap2 --version`;
	chomp $minimap_version;
}

my $samtools_version = `samtools --version | head -1 | sed 's/^samtools //'`;
chomp $samtools_version;

my $qualimap_version;
if(not defined(can_run("qualimap"))){
	$qualimap_version = "not detected";
}
else{
	$qualimap_version = `qualimap bamqc 2> /dev/null | head -4 | tail -1 | sed 's/^QualiMap v.//'`;
	chomp $qualimap_version;
}

my $bedtools_version;
if(not defined(can_run("bedtools"))){
	$bedtools_version = "not detected";
}
else{
	$bedtools_version = `bedtools --version | awk '{print \$2}' | sed 's/^v//'`;
	chomp $bedtools_version;
}

my $rscript_version;
if(not defined(can_run("Rscript"))){
	$rscript_version = "not detected";
}
else{
	$rscript_version = `Rscript --version 2>&1 | sed 's/^R scripting front-end version //;s/ .*//'`;
	chomp $rscript_version;
}

my $multiqc_version;
if(not defined(can_run("multiqc"))){
	$multiqc_version = "not detected";
}
else{
	$multiqc_version = `multiqc --version 2> /dev/null | awk '{print \$NF}'`;
	chomp $multiqc_version;
}

my $verbose_word = "No";
if($verbose == 1){
	$verbose_word = "Yes";
}

my $keep_tmp_word = "No";
if($keep_tmp == 1){
	$$keep_tmp_word = "Yes";
}

my $run_bamqc_switch_word = "Yes";
if($run_bamqc_switch == 0){
	$run_bamqc_switch_word = "No";
}

my $create_histo_switch_word = "Yes";
if($create_histo_switch == 0){
	$create_histo_switch_word = "No";
}

my $estimate_genome_size_switch_word = "Yes";
if($estimate_genome_size_switch == 0){
	$estimate_genome_size_switch_word = "No";
}

print "\n";
print "backmap.pl v$version\n";
print "\n";
print "Detected tools\n";
print "==============\n";
print "bwa:                  " . $bwa_version . "\n";
print "minimap2:             " . $minimap_version . "\n";
print "samtools:             " . $samtools_version . "\n";
print "qualimap:             " . $qualimap_version . "\n";
print "bedtools:             " . $bedtools_version . "\n";
print "Rscript:              " . $rscript_version . "\n";
print "multiqc:              " . $multiqc_version . "\n";
print "\n";
print "User defined input\n";
print "==================\n";
print "Output directory:     " . $out_dir . "\n";
if($assembly_path ne ""){
	print "Assembly:             " . $assembly_path . "\n";
	print "Paired reads:         ";
	print join("\n                      ",keys(%paired_filter)) . "\n";
	print "Unpaired reads:       ";
	print join("\n                      ",keys(%unpaired_filter)) . "\n";
	print "PacBio reads:         ";
	print join("\n                      ",keys(%pb_filter)) . "\n";
	print "Nanopore reads:       ";
	print join("\n                      ",keys(%ont_filter)) . "\n";
}
if(scalar(keys(%bam_filter)) > 0){
	print "Bam files:            ";
	print join("\n                      ",keys(%bam_filter)) . "\n";
}
print "Number of threads:    " . $threads . "\n";
print "Outpufile prefix:     " . $prefix . "\n";
print "Verbose:              " . $verbose_word . "\n";
print "Keep temporary files: " . $keep_tmp_word . "\n";
if($assembly_path ne ""){
	print "bwa options:          " . $bwa_opts . "\n";
	print "minimap2 options:     " . $minimap_opts . "\n";
}
print "Run qualimap bamqc:   " . $run_bamqc_switch_word . "\n";
if($run_bamqc_switch == 1){
	print "qualimap options:     " . $qm_opts . "\n";
}
print "Create cov histo:     " . $create_histo_switch_word . "\n";
print "Estimate genome size: " . $estimate_genome_size_switch_word . "\n";

if(defined $mkdir_cmd){
	exe_cmd($mkdir_cmd,$verbose,$dry);
}

if(defined $cmd){
	exe_cmd($cmd,$verbose,$dry);
}

#mapping will be executed if there are entries in %paired_filer or %unpaired_filter
#this happens if $assembly_path is set only
my $paired_counter = 0;
my @paired_bam = ();
foreach(keys(%paired_filter)){
	$paired_counter++;
	my ($for,$rev) = split(/,/,$_);
	$cmd = "bwa mem -t $threads $bwa_opts$index_path $for $rev 2> $out_dir/$prefix\_bwa_mem_paired$paired_counter.err | samtools view -1 -b - > $out_dir/$prefix.paired$paired_counter.bam";
	exe_cmd($cmd,$verbose,$dry);
	push(@paired_bam,"$out_dir/$prefix.paired$paired_counter.bam")
}

my $unpaired_counter = 0;
my @unpaired_bam = ();
foreach(keys(%unpaired_filter)){
	$unpaired_counter++;
	$cmd = "bwa mem -t $threads $bwa_opts$index_path $_ 2> $out_dir/$prefix\_bwa_mem_unpaired$unpaired_counter.err | samtools view -1 -b - > $out_dir/$prefix.unpaired$unpaired_counter.bam";
	exe_cmd($cmd,$verbose,$dry);
	push(@unpaired_bam,"$out_dir/$prefix.unpaired$unpaired_counter.bam");
}

my $pb_counter = 0;
my @pb_bam = ();
foreach(keys(%pb_filter)){
	$pb_counter++;
	$cmd = "minimap2 $minimap_opts-H -x map-pb -a -t $threads $assembly_path $_ 2> $out_dir/$prefix\_minimap_pb$pb_counter.err | samtools view -1 -b - > $out_dir/$prefix.pb$pb_counter.bam";
	exe_cmd($cmd,$verbose,$dry);
	push(@pb_bam,"$out_dir/$prefix.pb$pb_counter.bam");
}

my $ont_counter = 0;
my @ont_bam = ();
foreach(keys(%ont_filter)){
	$ont_counter++;
	$cmd = "minimap2 $minimap_opts-x map-ont -a -t $threads $assembly_path $_ 2> $out_dir/$prefix\_minimap_ont$ont_counter.err| samtools view -1 -b - > $out_dir/$prefix.ont$ont_counter.bam";
	exe_cmd($cmd,$verbose,$dry);
	push(@ont_bam,"$out_dir/$prefix.ont$ont_counter.bam");
}

my @merged_bam_file = ();
if($assembly_path ne ""){
	my $ill_bam_count = scalar(@paired_bam) + scalar(@unpaired_bam);
	
	my $paired_bam_files = join(" ",@paired_bam);
	my $unpaired_bam_files = join(" ",@unpaired_bam);
	my $pb_bam_files = join(" ",@pb_bam);
	my $ont_bam_files = join(" ",@ont_bam);
	
	if($ill_bam_count > 0){
		if($ill_bam_count == 1){
			my $single_bam = join(" ",@paired_bam,@unpaired_bam);
			$single_bam =~ s/^\s+//;
				$single_bam =~ s/\s+$//;
			$cmd = "ln -s $single_bam $out_dir/$prefix.bam";
			exe_cmd($cmd,$verbose,$dry);
			push(@merged_bam_file, "$out_dir/$prefix.bam");
		}
		else{
			$cmd = "samtools merge -@ $samtools_threads $out_dir/$prefix.bam $paired_bam_files $unpaired_bam_files";
			exe_cmd($cmd,$verbose,$dry);
			push(@merged_bam_file, "$out_dir/$prefix.bam");
		}
	}
	
	if(scalar(@pb_bam) > 0){
		if(scalar(@pb_bam) == 1){
			my $single_bam = $pb_bam[0];
			$cmd = "ln -s $single_bam $out_dir/$prefix.pb.bam";
			exe_cmd($cmd,$verbose,$dry);
			push(@merged_bam_file, "$out_dir/$prefix.pb.bam");
		}
		else{
			$cmd = "samtools merge -@ $samtools_threads $out_dir/$prefix.pb.bam $pb_bam_files";
			exe_cmd($cmd,$verbose,$dry);
			push(@merged_bam_file, "$out_dir/$prefix.pb.bam");
		}
	}
	
	if(scalar(@ont_bam) > 0){
		if(scalar(@ont_bam) == 1){
			my $single_bam = $ont_bam[0];
			$cmd = "ln -s $single_bam $out_dir/$prefix.ont.bam";
			exe_cmd($cmd,$verbose,$dry);
			push(@merged_bam_file, "$out_dir/$prefix.ont.bam");
		}
		else{
			$cmd = "samtools merge -@ $samtools_threads $out_dir/$prefix.ont.bam $ont_bam_files";
			exe_cmd($cmd,$verbose,$dry);
			push(@merged_bam_file, "$out_dir/$prefix.ont.bam");
		}
	}
	
}

my @sorted_bams;
if(scalar(keys(%bam_filter)) > 0){
	if($sort_bam_switch == 0){
		@sorted_bams = keys(%bam_filter);
	}
	if($sort_bam_switch == 1){
		@merged_bam_file = keys(%bam_filter);
	}
}

if($assembly_path ne "" or $sort_bam_switch == 1){
	foreach(@merged_bam_file){
		my $sorted_bam_file = $_;
		$sorted_bam_file =~ s/\.bam$/\.sort\.bam/;
		$cmd = "samtools sort -l 9 -@ $samtools_threads -T $out_dir/$prefix -o $sorted_bam_file $_";
		exe_cmd($cmd,$verbose,$dry);
		push(@sorted_bams,$sorted_bam_file);
	}
	
	if($keep_tmp == 0){
		my $tmp_bams = join(" ",@paired_bam,@unpaired_bam,@pb_bam,@ont_bam,@merged_bam_file);
		$cmd = "rm $tmp_bams";
		exe_cmd($cmd,$verbose,$dry);
	}
}

if($run_bamqc_switch == 1){
	if(defined(can_run("qualimap"))){
		foreach(@sorted_bams){
			my $bamqc_out = $_;
			$bamqc_out = (split(/\//,$bamqc_out))[-1];
			$bamqc_out =~ s/\.bam$/_stats/;
			$cmd = "qualimap bamqc $qm_opts-bam $_ -nt $threads -outdir $out_dir/$bamqc_out > $out_dir/$bamqc_out\_bamqc.log 2> $out_dir/$bamqc_out\_bamqc.err";
			exe_cmd($cmd,$verbose,$dry);
		}
		if(defined(can_run("multiqc")) and scalar(@sorted_bams) > 1){
			$cmd = "multiqc -s -m qualimap -o $out_dir $out_dir > $out_dir/multiqc.log 2> $out_dir/multiqc.err";
			exe_cmd($cmd,$verbose,$dry);
		}
	}
}

my %cov_files;
my %peak_cov;
my %n0_all;
my @global_ymax = ();

my $maxProcs = scalar(@sorted_bams);
if($threads < $maxProcs){
	$maxProcs = $threads;
}

if($create_histo_switch == 1){
	
	my $multiple_histos = Parallel::Loops->new($maxProcs);
	$multiple_histos->share(\%cov_files);
	$multiple_histos->share(\%peak_cov);
	$multiple_histos->share(\%n0_all);
	$multiple_histos->share(\@global_ymax);
	
	$multiple_histos->foreach( \@sorted_bams,  sub {
		
		my $tech = "Illumina";
		if($_ =~ m/\.pb\.sort\.bam$/){
			$tech = "PacBio";
		}
		if($_ =~ m/\.ont\.sort\.bam$/){
			$tech = "Nanopore";
		}
		
		my $filename = (split(/\//,$_))[-1];
		my $cov_hist_file = "$out_dir/$filename.cov-hist";
		$cmd = "samtools view -b -h -F 256 $_ | bedtools genomecov -ibam stdin -d | awk \'{print \$3}\' | sort -g | uniq -c | awk '{print \$2\"\\t\"\$1}' > $cov_hist_file";
		exe_cmd($cmd,$verbose,$dry);
		
		$cov_files{$tech} = $cov_hist_file;
		
		if($dry == 0){
			my $peak = `sort -rgk2 $cov_hist_file | awk \'\$1!=0{print \$1}\' | head -1`;
			chomp $peak;
			$peak_cov{$tech} = $peak;
			my $n0 = `awk \'\$1==0{print \$2}\' $cov_hist_file`;
			chomp $n0;
			$n0_all{$tech} = $n0;
			my $ymax = `sort -rgk2 $cov_hist_file | awk \'\$1!=0{print \$2}\' | head -1`;
			chomp $ymax;
			push(@global_ymax,$ymax);
			
			open(R,'>',$cov_hist_file . ".plot.r") or die "ERROR\tCould not open file " . $cov_hist_file . "plot.r\n";
			
			print R "x=read.table(\"$cov_hist_file\")\n";
			print R "pdf(\"$cov_hist_file.pdf\")\n";
			if($n0 < $ymax){
				print R "plot(x[,1],x[,2],log=\"x\",type=\"l\",xlab=\"Coverage\",ylab=\"Count\",main=\"$assembly\\n$tech\")\n";
			}
			else{
				print R "plot(x[,1],x[,2],ylim=c(0,$ymax),log=\"x\",type=\"l\",xlab=\"Coverage\",ylab=\"Count\",main=\"$assembly\\n$tech\")\n";
			}
			print R "text(2.5,$ymax,\"N(0)=$n0\")\n";
			print R "dev.off()\n";
			
			close R;
		}
		else{
			print STDERR "CMD\tsort -rgk2 $cov_hist_file | awk \'\$1!=0{print \$1}\' | head -1\n";
			print STDERR "CMD\tawk \'\$1==0{print \$2}\' $cov_hist_file\n";
			print STDERR "I would create $cov_hist_file.plot.r\n";
		}
		
		if(defined(can_run("Rscript"))){
			$cmd = "Rscript $cov_hist_file.plot.r > /dev/null 2> /dev/null";
			exe_cmd($cmd,$verbose,$dry);
		}
	});
	
	my @global_ymax = sort( {$b <=> $a} @global_ymax);
	
	if(scalar(@sorted_bams) > 1){
		my $rscript = "$out_dir/plot.all.r";
		if($prefix ne ""){
			$rscript = "$out_dir/$prefix.plot.all.r";
		}
		if($dry == 0){
			my @techs = ("Illumina","PacBio","Nanopore");
			
			open(RALL,'>',"$rscript") or die "ERROR\tCould not open file $rscript\n";
			
			for(my $i = 0; $i < scalar(@techs); $i++){
				if(exists($cov_files{$techs[$i]})){
					print RALL "$techs[$i]=read.table(\"$cov_files{$techs[$i]}\")\n";
				}
			}
			my $pdf = $rscript;
			$pdf =~ s/r$/pdf/;
			print RALL "pdf(\"$pdf\")\n";
			my @legend = ();
			my @lty = ();
			my @col = ();
			for(my $i = 0; $i < scalar(@techs); $i++){
				if(exists $cov_files{$techs[$i]}){
					push(@legend,"\"$techs[$i] N(0)=$n0_all{$techs[$i]}\"");
					push(@lty,"1");
					if($i == 0 and exists($cov_files{$techs[$i]})){
						print RALL "plot($techs[$i]\[,1],$techs[$i]\[,2],log=\"x\",type=\"l\",xlab=\"Coverage\",ylab=\"Count\",main=\"$assembly\",ylim=c(0,$global_ymax[0]))\n";
						push(@col,"\"black\"");
					}
					if($i == 1 and exists($cov_files{$techs[$i]})){
						print RALL "lines($techs[$i]\[,1],$techs[$i]\[,2],type=\"l\",col=\"blue\")\n";
						push(@col,"\"blue\"");
					}
					if($i == 2 and exists($cov_files{$techs[$i]})){
						print RALL "lines($techs[$i]\[,1],$techs[$i]\[,2],type=\"l\",col=\"red\")\n";
						push(@col,"\"red\"");
					}
				}
			}
			print RALL "legend(\"topright\",legend=c(" . join(",",@legend) . "),lty=c(" . join(",",@lty) . "),col=c(" . join(",",@col) . "))\n";
			print RALL "dev.off()\n";
			
			close RALL;
		}
		else{
			print STDERR "I would create $rscript\n";
		}
		
		if(defined(can_run("Rscript"))){
			$cmd = "Rscript $rscript > /dev/null 2> /dev/null";
			exe_cmd($cmd,$verbose,$dry);
		}
	}
	
}

if($estimate_genome_size_switch == 1){
	
	my %results;
		
        my $multiple_genome_size = Parallel::Loops->new($maxProcs);
        $multiple_genome_size->share(\%results);
	
	$multiple_genome_size->foreach( \@sorted_bams,  sub {
		if($dry == 1){
			print STDERR "CMD\tsort -rgk2 $_.cov-hist | awk \'\$1!=0{print \$1}\' | head -1\n";
		}
		
		$cmd = "samtools stats $_ > $_.stats 2> $_.stats.err";
		exe_cmd($cmd,$verbose,$dry);
		
		if($dry == 0){
			
			my $tech = "Illumina";
			if($_ =~ m/\.pb\.sort\.bam$/){
				$tech = "PacBio";
			}
			if($_ =~ m/\.ont\.sort\.bam$/){
				$tech = "Nanopore";
			}
			
			my $assembly_length = 0;
			my $assemly_perc = 0;
			if($assembly ne ""){
				open(IN,'<',"$assembly_path") or die "ERROR\tCould not open file " . $assembly_path . "\n";
				while(my $line = <IN>){
					chomp $line;
					if($line !~ m/^>/){
						$assembly_length = $assembly_length + length($line);
					}
				}
			}
			my $total_nucs = `grep "bases mapped (cigar):" $_.stats | awk -F'\\t' '{print \$3}'`;
			chomp $total_nucs;
			
			my $genome_size_estimate = $total_nucs / $peak_cov{$tech};

			$results{$tech} = $tech . " (" . $_ . ")\n" . "Mapped nucleotides:   " . round_format_pref($total_nucs) . "b\n" . "Peak coverage:        " . $peak_cov{$tech} . "\n" . "Genome size estimate: " . round_format_pref($genome_size_estimate) . "b\n";
			if($assembly_length > 0){
				$assemly_perc = sprintf("%.2f", ($assembly_length / $genome_size_estimate) * 100);
				$results{$tech} = $results{$tech} . "Assembly length:      " . round_format_pref($assembly_length) . "b ($assemly_perc% of estimate)\n";
			}
			
		}
		else{
			print STDERR "CMD\tgrep \"bases mapped (cigar):\" $_.stats | awk -F\'\\t\' \'{print \$3}\'\n";
		}
	});
	
	if($dry == 0){
		print "\n";
		print "Output\n";
		print "======\n";
		
		my @techs = ("Illumina","PacBio","Nanopore");
		for (my $i = 0; $i < scalar(@techs); $i++){
			if(exists($results{$techs[$i]})){
				print $results{$techs[$i]};
			}
		}
	}
}

exit;
