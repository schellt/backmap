#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';
use IPC::Cmd qw[can_run run];
use Number::FormatEng qw(:all);

my $version = "0.1";

sub print_help{
	print STDOUT "\n";
	print STDOUT "backmap.pl v$version\n";
	print STDOUT "\n";
	print STDOUT "Description:\n";
	print STDOUT "\tAutomatic mapping of paired and unpaired reads to an assembly, execution of\n\tqualimap bamqc and estimation of genome size from mapped nucleotides and\n\tpeak coverage.\n\tThe tools bwa, samtools, qualimap and bedtools need to be in your \$PATH.\n";
	print STDOUT "\n";
	print STDOUT "Usage:\n";
	print STDOUT "\tbackmap.pl [-a <assembly.fa> {-p <paired_1.fq>,<paired_2.fq> | -u <unpaired.fq>} |\n";
	print STDOUT "\t            -b <mapping.bam>]\n";
	print STDOUT "\n";
	print STDOUT "Mandatory:\n";
	print STDOUT "\t-a STR\t\tAssembly were reads should mapped to in fasta format\n";
	print STDOUT "\t-p STR\t\tTwo files with paired reads comma sperated\n";
	print STDOUT "\t\t\tCan be specified multiple times\n";
	print STDOUT "\t-u STR\t\tOne file with unpaired reads\n";
	print STDOUT "\t\t\tCan be specified multiple times\n";
	print STDOUT "\tOR\n";
	print STDOUT "\t-b STR\t\tBam file to calculate coverage from\n";
	print STDOUT "\t\t\tSkips read mapping\n";
	print STDOUT "\t\t\tOverrides -nh\n";
	print STDOUT "\n";
	print STDOUT "Options: [default]\n";
	print STDOUT "\t-o STR\t\tOutput directory [.]\n";
	print STDOUT "\t\t\tWill be created if not existing\n";
	print STDOUT "\t-t INT\t\tNumber of parallel executed processes [1]\n";
	print STDOUT "\t\t\tAffects bwa mem, samtools sort, qualimap bamqc\n";
	print STDOUT "\t-pre STR\tPrefix of output files [filename of -a or -b]\n";
	print STDOUT "\t-sort\t\tSort the bam file (-b) [off]\n";
	print STDOUT "\t-nq\t\tDo not run qualimap bamqc [off]\n";
	print STDOUT "\t-nh\t\tDo not create coverage histogram [off]\n";
	print STDOUT "\t\t\tImplies -ne\n";
	print STDOUT "\t-ne\t\tDo not estimate genome size [off]\n";
	print STDOUT "\t-kt\t\tKeep temporary bam files [off]\n";
	print STDOUT "\t-bo STR\t\tOptions passed to bwa [-a -c 10000]\n";
	print STDOUT "\t\t\tPass options with -bo \"<options>\"\n";
	print STDOUT "\t-qo STR\t\tOptions passed to qualimap [none]\n";
	print STDOUT "\t\t\tPass options with -qo \"<options>\"\n";
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

my $out_dir = abs_path("./");
my $assembly_path = "";
my $assembly = "";
#my $assembly_link = "";
my @paired = ();
my @unpaired = ();
my $threads = 1;
my $prefix = "";
my $verbose = 0;
my $bwa_opts = "-a -c 10000 ";
my $qm_opts = "";
my $create_histo_switch = 1;
my $estimate_genome_size_switch = 1;
my $run_bamqc_switch = 1;
my $keep_tmp = 0;
my $dry = 0;
my $bam = "";
my $sort_bam_switch = 0;
my $cmd;

my $input_error = 0;

if(scalar(@ARGV==0)){
	print_help;
}

for (my $i = 0; $i < scalar(@ARGV);$i++){
	if ($ARGV[$i] eq "-o"){
		$out_dir = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-a"){
		$assembly = (split /\//,$ARGV[$i+1])[-1];
		$assembly_path = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-p"){
		push(@paired,$ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-u"){
		push(@unpaired,$ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-b"){
		$bam = abs_path($ARGV[$i+1]);
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

if($assembly_path ne "" and $bam ne ""){
	print STDERR "ERROR\tSpecify either -a or -b\n";
	$input_error = 1;
}

if($assembly_path ne "" and $bam eq ""){
	if(not defined(can_run("bwa"))){
		print STDERR "ERROR\tbwa is not in your \$PATH\n";
		$input_error = 1;
	}
	if(not defined(can_run("samtools"))){
		print STDERR "ERROR\tsamtools is not in your \$PATH\n";
		$input_error = 1;
	}
}

if(not defined(can_run("qualimap"))){
	print STDERR "INFO\tqualimap is not in your \$PATH and will not be executed\n";
}

if(not defined(can_run("bedtools"))){
	print STDERR "INFO\tbedtools is not in your \$PATH\n";
	if($estimate_genome_size_switch == 1){
		print STDERR "WARNING\tGenome size estimation not possible\n";
		$estimate_genome_size_switch = 0;
	}
	if($bam ne ""){
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
	if(scalar(@paired) == 0 and scalar(@unpaired) == 0){
		print STDERR "ERROR\tNo reads specified!\n";
		$input_error = 1;
	}
}

if($bam ne ""){
	if(not -f $bam){
		print STDERR "ERROR\tFile $bam does not exist!\n";
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

if(not -d "$out_dir"){
	print STDERR "INFO\tCreating output directory $out_dir\n";
	$cmd="mkdir -p $out_dir";
	exe_cmd($cmd,$verbose,$dry);
}

if($prefix eq ""){
	if($assembly ne ""){
		$prefix = $assembly;
	}
	if($bam ne ""){
		$prefix = (split /\//,$bam)[-1];
	}
	print STDERR "INFO\tSetting prefix to $prefix\n";
}

if($bam ne "" and defined(can_run("bedtools"))){
	if($create_histo_switch == 0 or $estimate_genome_size_switch == 0){
		$create_histo_switch = 1;
		print STDERR "INFO\tCreating coverage histogram\n";
	}
}

my %paired_filter;

if($assembly_path ne ""){
	foreach(@paired){
		my @pair = split(/,/,$_);
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

if($assembly_path ne ""){
	if(scalar(keys(%paired_filter)) == 0 and scalar(keys(%unpaired_filter)) == 0){
		print STDERR "ERROR\tNo existing read files specified!\n";
		exit 1;
	}
}

my $index_path = $out_dir . "/" . $assembly;

if($assembly_path ne ""){
	if(not -f $assembly_path . ".amb" or not -f $assembly_path . ".ann" or not -f $assembly_path . ".bwt" or not -f $assembly_path . ".pac" or not -f $assembly_path . ".sa"){
		$cmd = "bwa index -p $index_path $assembly_path > $out_dir/$prefix\_bwa_index.log 2> $out_dir/$prefix\_bwa_index.err";
	}
	else{
		print STDERR "INFO\tIndex files already existing\n";
		$index_path = $assembly_path;
	}
}

my $bwa_version = `bwa 2>&1 | head -3 | tail -1 | sed 's/^Version: //'`;
chomp $bwa_version;

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
print "Detected tools\n";
print "==============\n";
print "bwa:                  " . $bwa_version . "\n";
print "samtools:             " . $samtools_version . "\n";
print "qualimap:             " . $qualimap_version . "\n";
print "bedtools:             " . $bedtools_version . "\n";
print "Rscript:              " . $rscript_version . "\n";
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
}
if($bam ne ""){
	print "Bam file:             " . $bam . "\n";
}
print "Number of threads:    " . $threads . "\n";
print "Outpufile prefix:     " . $prefix . "\n";
print "Verbose:              " . $verbose_word . "\n";
print "Keep temporary files: " . $keep_tmp_word . "\n";
if($assembly_path ne ""){
	print "bwa options:          " . $bwa_opts . "\n";
}
print "Run qualimap bamqc:   " . $run_bamqc_switch_word . "\n";
if($run_bamqc_switch == 1){
	print "qualimap options:     " . $qm_opts . "\n";
}
print "Create cov histo:     " . $create_histo_switch_word . "\n";
print "Estimate genome size: " . $estimate_genome_size_switch_word . "\n";

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

my $merged_bam_file;
if($assembly_path ne ""){
	my $bam_count = scalar(@paired_bam) + scalar(@unpaired_bam);
	
	my $paired_bam_files = join(" ",@paired_bam);
	my $unpaired_bam_files = join(" ",@unpaired_bam);
	$merged_bam_file = $out_dir . "/" . $prefix . ".bam";
	
	if($bam_count == 1){
		my $single_bam = join(" ",@paired_bam,@unpaired_bam);
		$single_bam =~ s/^\s+//;
		$single_bam =~ s/\s+$//;
		$cmd = "ln -s $single_bam $merged_bam_file";
		exe_cmd($cmd,$verbose,$dry);
	}
	else{
		$cmd = "samtools merge $merged_bam_file $paired_bam_files $unpaired_bam_files";
		exe_cmd($cmd,$verbose,$dry);
	}
}

my $sorted_bam_file;
if($bam ne ""){
	if($sort_bam_switch == 0){
		$sorted_bam_file = $bam;
	}
	if($sort_bam_switch == 1){
		$merged_bam_file = $bam;
		$sorted_bam_file = $out_dir . "/" . $prefix . ".sort.bam";
	}
}

if($assembly_path ne ""){
	$sorted_bam_file = $out_dir . "/" . $prefix . ".sort.bam";
}

if($assembly_path ne "" or $sort_bam_switch == 1){
	$cmd = "samtools sort -l 9 -@ $threads -T $prefix -o $sorted_bam_file $merged_bam_file";
	exe_cmd($cmd,$verbose,$dry);
	
	if($keep_tmp == 0){
		my $tmp_bams = join(" ",@paired_bam,@unpaired_bam);
		$cmd = "rm $merged_bam_file $tmp_bams";
		exe_cmd($cmd,$verbose,$dry);
	}
}

if($run_bamqc_switch == 1){
	if(defined(can_run("qualimap"))){
		$cmd = "qualimap bamqc $qm_opts-bam $sorted_bam_file -nt $threads > $out_dir/$prefix\_bamqc.log 2> $out_dir/$prefix\_bamqc.err";
		exe_cmd($cmd,$verbose,$dry);
	}
}

my $cov_hist_file;
my $peak_cov;
if($create_histo_switch == 1){
	$cov_hist_file = "$out_dir/$prefix.cov-hist";
	
	$cmd = "samtools view -b -h -F 256 $sorted_bam_file | bedtools genomecov -ibam stdin -d | awk \'{print \$3}\' | sort -g | uniq -c | awk '{print \$2\"\\t\"\$1}' > $cov_hist_file";
	exe_cmd($cmd,$verbose,$dry);

	if($dry == 0){
		$peak_cov = `sort -rgk2 $cov_hist_file | awk \'\$1!=0{print \$1}\' | head -1`;
		chomp $peak_cov;
		my $n0 = `awk \'\$1==0{print \$2}\' $cov_hist_file`;
		chomp $n0;
		my $ymax = `sort -rgk2 $cov_hist_file | awk \'\$1!=0{print \$2}\' | head -1`;
		chomp $ymax;
	
		open(R,'>',$cov_hist_file . ".plot.r") or die "ERROR\tCould not open file " . $cov_hist_file . "plot.r\n";
		
		print R "x=read.table(\"$cov_hist_file\")\n";
		print R "pdf(\"$cov_hist_file.pdf\")\n";
		if($n0 < $ymax){
			print R "plot(x[,1],x[,2],log=\"x\",type=\"l\",xlab=\"Coverage\",ylab=\"Count\")\n";
		}
		else{
			print R "plot(x[,1],x[,2],ylim=c(0,$ymax),log=\"x\",type=\"l\",xlab=\"Coverage\",ylab=\"Count\")\n";
		}
		print R "text(2.5,$ymax,\"N(0)=$n0\")\n";
		print R "dev.off()";
		
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
}

if($estimate_genome_size_switch == 1){
	
	if($dry == 0){
		$peak_cov = `sort -rgk2 $cov_hist_file | awk \'\$1!=0{print \$1}\' | head -1`;
		chomp $peak_cov;
	}
	else{
		print STDERR "CMD\tsort -rgk2 $cov_hist_file | awk \'\$1!=0{print \$1}\' | head -1\n";
	}

	$cmd = "samtools stats $sorted_bam_file > $out_dir/$prefix.stats 2> $out_dir/$prefix.stats.err";
	exe_cmd($cmd,$verbose,$dry);
	
	my $total_nucs;	
	if($dry == 0){
		$total_nucs = `grep "bases mapped (cigar):" $out_dir/$prefix.stats | awk -F'\\t' '{print \$3}'`;
	}
	else{
		print STDERR "CMD\tgrep \"bases mapped (cigar):\" $out_dir/$prefix.stats | awk -F\'\\t\' \'{print \$3}\'\n";
	}
	
	if($dry == 0){
		print "\n";
		print "Output\n";
		print "==================\n";
		print "Mapped nucleotides:   " . format_pref($total_nucs) . "b\n";
		print "Peak coverage:        " . $peak_cov . "\n";
		
		my $genome_size_estimate = $total_nucs / $peak_cov;
		print "Genome size estimate: " . format_pref($genome_size_estimate) . "b\n";
	}
}

exit;
