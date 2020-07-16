#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path;

my $assembly_type = "";
my $single_file = "";
my $forward_file = "";
my $reverse_file = "";
my $interleaved_file = "";
my $param_file = "";
my $out_dir = "";
my $command;
my $program = $0;
my $dir = substr($0, 0, length($0)-11);
my $start = time();
my $starttime;
my $endtime;
my $help;
my $MOL=40;
my $OL=40;
my $assembly_file = "";
my $read_tag = "";
my $mergedFastq = "";
my $ncpu ;
my $max_ext_len = 150;
my $fgs_complete;
my $fgs_train_reads ;
my $fgs_train_other ;
my $fgs_drop_indel ;
my $mga_genome ;
my $sga_preprocess ;
my $sga_index ;
my $sga_index_algorithm ;
my $sga_overlap;
my $sga_overlap_min_len;
my $sga_rmdup ;
my $sga_assemble ;
my $sga_assemble_min_len ;
my $sga_merge;
my $sga_bwt2fa ;
my $sga_correct ;
my $sga_correct_algo ;
my $sga_correct_kmer;
my $sga_correct_min_olap ;
my $sga_fmmerge ;
my $sga_fmmerge_min_olap ;
my $sga_filter ;
my $sga_filter_kmer ;
my $sga_filter_kmer_thresh ;
my $spades;
my $spades_only_assembler;
my $spades_error_correct;
my $spades_mem_limit ;
my $spades_disable_gzip;
my $spades_continue ;
my $spades_kmer;
my $spades_rr;
my $bwa ;
my $bwa_min_seed ;
my $bwa_min_score;
my $plass;
my $plass_assemble;
my $plass_num_iter;
my $plass_min_len;
my $FILEA;
my $FILEB;
my $val;

#=====================================================================
### Print program usage
##====================================================================
sub print_usage{

    print "\nUSAGE: ./impp_run.pl [options] -o <output_directory>";
    print"\n\nInput data options:\n\n";
    print " --12            <filename>    : fastq file with interlaced forward and reverse paired-end reads\n";
    print " -1              <filename>    : fastq file with forward paired-end reads\n";
    print " -2	         <filename>    : fastq file with reverse paired-end reads\n";
    print " -s	         <filename>    : fastq file with unpaired reads\n";
    print " -a/--assembler  <0 or 1>      : [optional] type of assembler to run (0: sga, 1: spades)\n";
    print " -o/--outdir     <dirname>     : [required] output directory name including the full path\n";
    print " -m/--max-len    <int>         : [optional] maximum extension length for anchors (default: 150) \n";
    print " -p/--param-file <filename>    : [optional] parameter file\n";
    print " -h/--help                     : print help message\n";
    print "\n\n";
    print " NOTE: Input FASTQ read file must be specified in either --12, -1 & -2 or -s format. \n\n"
}


# Getting program options
my($reads_fa, $reads_fq);
getOpts();

#Calling FGS on reads
callFGS_Reads();

#Calling SGA/SPAdes on reads
callAssembly();

# Running impp to find and extend anchors
callIMPP();

# Calling FGS on pathset
callFGS_Paths();

# Running bwa on edges and paths
callBWA();

# Generating 6-frame translation on the uncalled reads from pathset
callSixFrameFilter();

# Calling peptide assembler (Plass)
callPeptideAssembler();

# Clean up and generate final output of iMPP
cleanOutput();

$endtime = time();
print "Successfully completed running iMPP in : ";
getElapsedTime($endtime - $start);
print "\n";

#====================================================================
## Process the input parameter options
#====================================================================
sub getOpts
{
GetOptions(
          's=s'             =>    \$single_file,
          '1=s'             =>    \$forward_file,
          '2=s'             =>    \$reverse_file,
          '12=s'            =>    \$interleaved_file,
          'assembly|a=s'    =>    \$assembly_type,
          'outdir|o=s'      =>    \$out_dir,
          'max-len|m=s'     =>    \$max_ext_len,
          'param-file|p=s'  =>    \$param_file,
          'help|h'          =>    \$help
           );

if ( $help ) {
    print_usage();
     exit;
}

if (length($single_file)==0 && length($forward_file)==0 && length($reverse_file)==0 && length($interleaved_file)==0){
    print "\nERROR: Input FASTQ read(s) file was not specified.\n\n";
    print_usage();
    exit;
}
else{
  if (! -e $single_file){
    if(! -e $interleaved_file){
      if( (! -e $forward_file) || (! -e $reverse_file) ){
        print "\nERROR: Both forward and reverse reads must be specified. See --help for details.\n\n";
      }
      else{
	my($FILEA, $FILEB);
        mergeFastqFiles();
      }
    }
    else{
      #$read_tag = substr($interleaved_file, 0, rindex ($interleaved_file, '.'));
      $command = "python ".$dir."utils/fastq2fasta.py $interleaved_file $out_dir";
      system($command);
    }
  }
  else{
    #$read_tag = substr($single_file, 0, rindex ($single_file, '.'));
    $command = "python ".$dir."utils/fastq2fasta.py $single_file $out_dir";
    system($command);

  }
}
$reads_fa = $out_dir."/reads.fq2fa.fasta";
$reads_fq = $out_dir."/reads.re.fastq";

if (length($assembly_type)==0){ $assembly_type = "0";	}
else{
    unless ($assembly_type eq "0" || $assembly_type eq "1"){
    	print "\nERROR: Please enter a valid input assembler. 0 for sga, 1 for spades. \n\n";
    	print_usage();
    	exit;
	}
}

if (length($out_dir) == 0 ){
    print "\nERROR: An output directory name must be specified.\n\n";
    print_usage();
    exit;
}


if ($max_ext_len && $max_ext_len >= 300 ){
    print "\nERROR: Maximum extend length threshold too high. Please select a lower value for optimal results. [default=300] \n\n";
    print_usage();
    exit;
}

if (length($param_file) == 0 ){
  if(! -e $dir."params/parameters.txt" ){
    print "\nERROR: Parameter file missing. Make sure parameters.txt file exists in $dir/params/ directory.\n";
  }
  else {  $param_file = $dir."params/parameters.txt"; }
}

#Loading parameters from parameter file
#loadParameters();

}


#====================================================================
# Call ORFs on input reads using FGS
#====================================================================
sub callFGS_Reads
{
  my $starttime = time();
  print STDERR "\n-----------------------------------------\n";
  print STDERR "\nRunning FragGeneScan on input reads ...\n";
  print STDERR "\n-----------------------------------------\n";
  $command = $dir."lib/FragGeneScan1.31/run_FragGeneScan.pl";
  $command .= " -genome=$reads_fa";
  $command .= " -out=$out_dir/reads";
  $command .= " -complete=0" ;
  $command .= " -train=illumina_5";
  $command .= " -thread=16";
  print "$command\n";
  system($command);

  $endtime = time();
  print STDERR "\n-----------------------------------------\n";
  print "FragGeneScan run on reads completed in : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";

}

#====================================================================
# Call Assembler ( String graph or de Bruijn Graph)
#====================================================================
sub callAssembly
{
  $starttime = time();
  if($assembly_type == 0)
    {
      #print STDERR "\n-----------------------------------------\n";
      print STDERR "\nRunning SGA assembler on input reads ...\n";
      print STDERR "\n-----------------------------------------\n";
      ## call SGA assembler
      my $SGA_BIN= $dir."lib/sga";
      print $dir;
      $command = "$SGA_BIN preprocess -o assembly.fq $reads_fa";
      system($command);
      $command= "$SGA_BIN index -a ropebwt -t 16 assembly.fq";
      system($command);
      $command = "$SGA_BIN rmdup -t 16 assembly.fq";
      system($command);
      $command = "$SGA_BIN overlap -m $MOL assembly.rmdup.fa";
      system($command);
      #$command = "$SGA_BIN assemble -m $OL -g 0 --transitive-reduction -r 10 -o assembly.m$OL $out_dir/assembly.rmdup.asqg.gz";
      #system($command);
      $command = "mv assembly.* $out_dir";
      system($command);


      ## Generating String Graph from the overlap graph output by SGA
      my $olap_file = $out_dir."/assembly.rmdup.asqg.gz";
      print $olap_file;
      if (-e $olap_file) {
          $command = "gunzip $olap_file";
          system($command);
        }
        else {
          print "\n Cannot open $olap_file";
        }
        $command = $dir."bin/impp_getSG";
        $command .= " -g $out_dir/assembly.rmdup.asqg";
        system($command);
        $assembly_file = $out_dir."/assembly.rmdup.StringGraph.fq";

        $endtime = time();

        print STDERR "\n-----------------------------------------\n";
        print STDERR "\nSGA assembly run completed in : ";
        getElapsedTime($endtime - $starttime);
        print STDERR "\n-----------------------------------------\n";

        # Calling FGS on assembly edges
        callFGS_Edges();
    }
    elsif ($assembly_type == 1)
    {
      ## call SPAdes assembler
      my $SPAdes_bin= $dir."lib/SPAdes-3.9.0-Linux/bin/metaspades.py";
      $command = "$SPAdes_bin --only-assembler -t 16 -m 500 --12 $reads_fa -o $out_dir/assembly"; 
      system($command);

      ## modifying the SPAdes graph
      $command = "python ";
      $command .= $dir."src/renameDB.py";
      $command .= " $out_dir/assembly/assembly_graph.fastg";
      system($command);
      $assembly_file = $out_dir."/assembly/assembly_graph.dbGraph.fq";

      $endtime = time();
      print STDERR "\n-----------------------------------------\n";
      print STDERR "\nSPAdes assembly run completed in : ";
      getElapsedTime($endtime - $starttime);
      print STDERR "\n-----------------------------------------\n";

      # Calling FGS on assembly edges
      callFGS_Edges();

    }
    else{
    print "ERROR: Input assembler value is invalid. Please select either 0 (for sga) or 1 (for spades) .\n"; }

}

#====================================================================
# Call ORFs on assembly graph edges using FGS
#====================================================================
sub callFGS_Edges
{
  $starttime = time();
 #print STDERR "\n-----------------------------------------\n";
  print STDERR "\nRunning FragGeneScan on edges ...\n";
  print STDERR "\n-----------------------------------------\n";
  $command = $dir."lib/FragGeneScan1.31/"."run_FragGeneScan.pl";
  $command .= " -genome=$assembly_file";
  $command .= " -out=$out_dir/edges";
  $command .= " -complete=0" ;
  $command .= " -train=complete";
  $command .= " -thread=16";
  print "$command\n";
  system($command);

  $endtime = time();
  print STDERR "\n-----------------------------------------\n";
  print "FragGeneScan run on edges completed in : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";

}


#=====================================================================
## Call ORFs on paths  using FGS
##====================================================================
sub callFGS_Paths
{
$starttime = time();
#print STDERR "\n-----------------------------------------\n";
print STDERR "\nRunning FragGeneScan on candidate paths generated by iMPP ...\n";
print STDERR "\n-----------------------------------------\n";
$command = $dir."lib/FragGeneScan1.31/"."run_FragGeneScan.pl";
$command .= " -genome=$out_dir/paths.fa";
$command .= " -out=$out_dir/paths";
$command .= " -complete=0" ;
$command .= " -train=complete";
$command .= " -thread=16";
print "$command\n";
system($command);

$endtime = time();
print STDERR "\n-----------------------------------------\n";
print "FragGeneScan run on paths completed in : ";
getElapsedTime($endtime - $starttime);
print STDERR "\n-----------------------------------------\n";
}

#======================================================================
### Aligning reads to edges and paths using BWA
###====================================================================
sub callBWA
{
my $starttime = time();
print STDERR "\nRunning BWA index (edges)...\n";
$command = $dir."lib/bwa-0.7.17/"."bwa index";
$command .= " -p $out_dir/edges.index";
$command .= " $assembly_file";
print "$command\n";
system($command);

print STDERR "\nRunning BWA mem...\n";
$command = $dir."lib/bwa-0.7.17/"."bwa mem";
$command .= " $out_dir/edges.index";
$command .= " $reads_fa";
$command .= " -T 90 -t 16";
$command .= " -o $out_dir/edges.sam";
print "$command\n";
system($command);
$endtime = time();
print "Time taken to complete bwa run on edges: \t";
getElapsedTime($endtime - $starttime);
print "\n";

$starttime = time();
print STDERR "\nRunning BWA index (paths)...\n";
$command = $dir."lib/bwa-0.7.17/"."bwa index";
$command .= " -p $out_dir/paths.index";
$command .= " $out_dir/paths.fa";
print "$command\n";
system($command);

print STDERR "\nRunning BWA mem...\n";
$command = $dir."lib/bwa-0.7.17/"."bwa mem";
$command .= " $out_dir/paths.index";
$command .= " $reads_fa";
$command .= " -T 90 -t 16";
$command .= " -o $out_dir/paths.sam";
print "$command\n";
system($command);

$endtime = time();
print "Time taken to complete bwa run on paths : ";
getElapsedTime($endtime - $starttime);
print "\n";
}


#=====================================================================
## Calling impp - finding anchors and extending Anchors to get paths
##====================================================================
sub callIMPP
{
  $starttime = time();
  print STDERR "\n Extending anchors..\n";

  my $command = $dir."bin/impp_extendAnchor";
  $command .= " -g $assembly_file ";
  $command .= " -a $out_dir/edges.gff";
  $command .= " -l $max_ext_len ";
  $command .= " -p $out_dir/paths.fa";

  system( "$command" ) == 0 or die "impp call failed: $?";

  $endtime = time();

  print "Time taken to find and extend Anchors : ";
  getElapsedTime($endtime - $starttime);
}


#=====================================================================
## Generating six frame translation
###====================================================================
sub callSixFrameFilter
{
  $starttime = time();
  #print STDERR "\n-----------------------------------------\n";
  print STDERR "\nRunning six frame translation filter on any un-called reads..\n";
  print STDERR "\n-----------------------------------------\n";
  my $command = "python ";
  $command .= $dir."utils/sixFrameFilter.py";
  $command .= " $out_dir/reads.gff ";
  $command .= " $out_dir/edges.sam";
  $command .= " $out_dir/paths.sam ";
  $command .= " $reads_fa";
  $command .= " $out_dir";
  system( "$command" ) == 0 or die "iMPP six frame reads filter failed: $?";
  $endtime = time();

  #print STDERR "\n-----------------------------------------\n";
  print "Time taken to generate 6 frame translation : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";
}


#=====================================================================
### Generating final peptide assembly
####====================================================================
sub callPeptideAssembler
{

  $starttime = time();
  print STDERR "\n Generating peptide assembly..\n";
  my $command = $dir."lib/plass/bin/plass assemble";
  $command .= " --threads 16 --min-length 30";
  $command .= " $out_dir/reads.filter.fasta $out_dir/assembled_proteins.faa";
  $command .= " $out_dir/plass.tmp ";

  system( "$command" ) == 0 or die "impp peptide assembly step failed: $?";
  system("rm -rf plass.tmp");

  $endtime = time();

  print STDERR "\n-----------------------------------------\n";
  print "Time taken to complete peptide assembly : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";
}

#==================================================
# Program to merge seperate forward and reverse paired end read files into single interleaved file
sub mergeFastqFiles
{
open ($FILEA,'<', $forward_file);
open ($FILEB, '<' ,$reverse_file);
open (my $mergedFastq, '>', $out_dir."/reads.merged.fq");
while(<$FILEA>) {
    print $mergedFastq $_;
    $_ = <$FILEA>;
    print $mergedFastq $_;
    $_ = <$FILEA>;
    print $mergedFastq $_;
    $_ = <$FILEA>;
    print $mergedFastq $_;

    $_ = <$FILEB>;
    print $mergedFastq $_;
    $_ = <$FILEB>;
    print $mergedFastq $_;
    $_ = <$FILEB>;
    print $mergedFastq $_;
    $_ = <$FILEB>;
    print $mergedFastq $_;
}

close($FILEA);
close($FILEB);
#print $mergedFastq;
#Program to convert fastq to fasta
$command = "python ".$dir."utils/fastq2fasta.py ".$out_dir."/reads.merged.fq";
print $command;
system($command);

}

#=====================================================================
#### Loading parameters for running 3rd party softwares
#####=================================================================
=begin comment
sub loadParameters
{
	open(IN, "< $param_file") or die $!;
	while ( <IN>) {
		chomp;
		next if /^#/;
		my ( $opt, $val ) = split;
		for ( $opt ) {
			when ("thread")                    {$ncpu     = $val; }
			when ("mlen")                      {$max_ext_len = $val;}
			when ("fgs_complete")              {$fgs_complete = $val; }
			when ("fgs_train_reads")           {$fgs_train_reads    = $val; }
      			when ("fgs_train_other")           {$fgs_train_other    = $val; }
			when ("fgs_dropindel")             {$fgs_drop_indel   = $val; }
			when ("assemble")                  {$mga_genome   = $val; }
      when ("sga_preprocess")             {$sga_preprocess = $val;}
      when ("sga_index")                  {$sga_index = $val;}
      when ("sga_index_algorithm ")       {$sga_index_algorithm = $val;}
      when ("sga_overlap")                {$sga_overlap = $val;}
      when ("sga_overlap_min_len")        {$sga_overlap_min_len = $val;}
      when ("sga_rmdup")                  {$sga_rmdup = $val;}
      when ("sga_assemble")               {$sga_assemble = $val;}
      when ("sga_assemble_min_len")       {$sga_assemble_min_len = $val;}
      when ("sga_merge")                  {$sga_merge = $val;}
      when ("sga_bwt2fa")                 {$sga_bwt2fa = $val;}
      when ("sga_correct")                {$sga_correct = $val;}
      when ("sga_correct_algo")           {$sga_correct_algo = $val;}
      when ("sga_correct_kmer")           {$sga_correct_kmer = $val;}
      when ("sga_correct_min_olap")       {$sga_correct_min_olap = $val;}
      when ("sga_fmmerge")                {$sga_fmmerge = $val;}
      when ("sga_fmmerge_min_olap")       {$sga_fmmerge_min_olap = $val;}
      when ("sga_filter")                 {$sga_filter = $val;}
      when ("sga_filter_kmer")            {$sga_filter_kmer = $val;}
      when ("sga_filter_kmer_thresh")     {$sga_filter_kmer_thresh = $val;}
      when ("spades")                     {$spades = $val;}
      when ("spades_only_assembler")      {$spades_only_assembler = $val;}
      when ("spades_error_correct")       {$spades_error_correct = $val;}
      when ("spades_memory_limit")        {$spades_mem_limit = $val;}
      when ("spades_disable_gzip")        {$spades_disable_gzip = $val;}
      when ("spades_continue")            {$spades_continue = $val;}
      when ("spades_kmer")                {$spades_kmer = $val;}
      when ("spades_rr")                  {$spades_rr = $val;}
      when ("bwa")                        {$bwa = $val;}
      when ("bwa_min_seed")               {$bwa_min_seed = $val;}
      when ("bwa_min_score")              {$bwa_min_score = $val;}
      when ("plass")                      {$plass = $val;}
      when ("plass_assemble")             {$plass_assemble = $val;}
      when ("plass_num_iter")             {$plass_num_iter = $val;}
      when ("plass_min_length")           {$plass_min_len = $val;}
		}
	}
	close(IN);
	$out_dir = File::Spec->rel2abs($out_dir);

  checkParameters();
}

#=====================================================================
#### Program to check input parameters in parameter file
#####=================================================================
sub checkParameters
{
	my $max_cpu = `grep -w "processor" /proc/cpuinfo | wc -l`; chomp $max_cpu;
	if ( $ncpu > $max_cpu ) {$ncpu = $max_cpu; }
	else { die "ERROR: Invalid value for number of threads option\n";}
	if ( ! defined $max_ext_len ) { die "ERROR: Max extend length must be given\n"; }
	if ( $fgs_complete != 0 && $fgs_complete != 1 ) { die "ERROR: Invalid value for FGS complete option\n"; }
	if ( $fgs_train_reads ne "complete"   &&
		 $fgs_train_reads ne "sanger_5"   &&
		 $fgs_train_reads ne "sanger_10"  &&
		 $fgs_train_reads ne "454_10"     &&
		 $fgs_train_reads ne "454_30"     &&
		 $fgs_train_reads ne "illumina_5" &&
		 $fgs_train_reads ne "illumina_10" ) { die "ERROR: Invalid value for FGS train option for input reads\n"; }
  if ( $fgs_train_other ne "complete"   &&
   		 $fgs_train_other ne "sanger_5"   &&
   		 $fgs_train_other ne "sanger_10"  &&
   		 $fgs_train_other ne "454_10"     &&
   		 $fgs_train_other ne "454_30"     &&
   		 $fgs_train_other ne "illumina_5" &&
   		 $fgs_train_other ne "illumina_10" ) { die "ERROR: Invalid value for FGS train option for processing\n"; }
	if ( $fgs_drop_indel != 0 && $fgs_drop_indel != 1 ) { die "ERROR: Invalid value for FGS dropindel option\n"; }
	if ( $assemble != 0 && $assemble != "1" ) { die "ERROR: Invalid value for assembler\n"; }
	if ( $sga_index != 0 && $sga_index != 1 ) { die "ERROR: Invalid sga index option\n"; }
  if ( $sga_index_algorithm != 0 && $sga_index_algorithm != 1 ) { die "ERROR: Invalid sga index algorithm option\n"; }
  if ( $sga_overlap != 0 && $sga_overlap != 1 ) { die "ERROR: Invalid sga overlap option\n"; }
  if ( $sga_overlap_min_len < 1 ) { die "ERROR: Invalid value sga overlap min length\n"; }
  if ( $sga_rmdup != 0 && $sga_rmdup != 1 ) { die "ERROR: Invalid sga rmdup option\n"; }
  if ( $sga_assemble!= 0 && $sga_assemble != 1 ) { die "ERROR: Invalid sga assemble option\n"; }
  if ( $sga_assemble_min_len < 1 ) { die "ERROR: Invalid value sga assemble min length\n"; }
	if ( $sga_merge != 0 && $sga_merge != 1 ) { die "ERROR: Invalid sga merge option\n"; }
  if ( $sga_bwt2fa != 0 && $sga_bwt2fa != 1 ) { die "ERROR: Invalid sga bwt2fa option\n"; }
  if ( $sga_correct != 0 && $sga_correct != 1 ) { die "ERROR: Invalid sga correct option\n"; }
  if ( $sga_correct_algo != 0 && $sga_correct_algo != 1 && $sga_correct_algo !=2 ) { die "ERROR: Invalid sga correct algorithm option\n"; }
  if ( $sga_correct_kmer < 1 ) { die "ERROR: Invalid value sga correct kmer length\n"; }
  if ( $sga_correct_min_olap < 1 ) { die "ERROR: Invalid value sga correct min overlap length\n"; }
  if ( $sga_fmmerge != 0 && $sga_fmmerge != 1 ) { die "ERROR: Invalid sga fmmerge option\n"; }
  if ( $sga_fmmerge_min_olap < 1 ) { die "ERROR: Invalid value sga fmmerge min overlap length\n"; }
  if ( $sga_filter != 0 && $sga_filter != 1 ) { die "ERROR: Invalid sga filter option\n"; }
  if ( $sga_filter_kmer < 1 ) { die "ERROR: Invalid value sga filter kmer length\n"; }
  if ( $sga_filter_kmer_thresh < 1 ) { die "ERROR: Invalid value sga filter kmew threshold\n"; }
  if ($spades ne "s") { die "ERROR: Invalid spades option\n";}
  if ( $spades_only_assembler != 0 && $spades_only_assembler != 1 ) { die "ERROR: Invalid spades assembler option\n"; }
  if ( $spades_error_correct != 0 && $spades_error_correct != 1 ) { die "ERROR: Invalid spades error correct option\n"; }
  if ( $spades_continue != 0 && $spades_continue != 1 ) { die "ERROR: Invalid spades continue option\n"; }
  if ( $spades_disable_gzip != 0 && $spades_disable_gzip != 1 ) { die "ERROR: Invalid spades disable gzip option\n"; }
  if ( $spades_rr != 0 && $spades_rr != 1 ) { die "ERROR: Invalid spades repeat resolution option\n"; }
  if ( $spades_mem_limit < 1 ) { die "ERROR: Invalid spades memory limit \n"; }
  if ( $bwa_min_seed < 0 ) { die "ERROR: Invalid bwa min seed option\n"; }
  if ( $bwa_min_score < 0  { die "ERROR: Invalid bwa min score option\n"; }
  if ( $plass != 0 && $plass != 1 ) { die "ERROR: Invalid plass option\n"; }
  if ( $plass_assemble != 0 && $plass_assemble != 1 ) { die "ERROR: Invalid plass assemble option\n"; }
  if ( $plass_min_len < 0 ) { die "ERROR: Invalid plass min length option\n"; }
  if ( $plass_num_iter < 0 ) { die "ERROR: Invalid plass num iterations option\n"; }


}
=end comment
=cut


#=====================================================================
### Clean up final output directory
###====================================================================
sub cleanOutput
{
  my $plass_dir = "$out_dir/plass.tmp/";
  unlink glob "'$out_dir/reads*.*'";  
  unlink glob "'$out_dir/edges*.*'"; 
  unlink glob "'$out_dir/assembly*.*'"; 
  rename("$out_dir/paths.gff", "$out_dir/orfs.gff");
  rename("$out_dir/paths.faa", "$out_dir/orfs.faa");
  rename("$out_dir/paths.ffn", "$out_dir/orfs.ffn");
  unlink glob "'$out_dir/paths*.*'"; 
  rmtree $plass_dir;
  
}
#=====================================================================
## Get total elapsed time
##====================================================================
sub getElapsedTime
{
	my $input = $_[0];
	my $hour;
	my $min;
	my $sec;
	my $str;

	$sec = $input % 60;
	$input = int($input / 60);
	$min = $input % 60;
	$input = int($input / 60);
	$hour = $input;

	$str = $hour . " hours " . $min . " minutes and " . $sec . " seconds.\n";
	print $str;
	print "\n";
}
