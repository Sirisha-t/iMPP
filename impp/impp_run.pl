#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path;

my @exts = qw(.fq .fastq);
my $assembly_type = "";
my $single_file = "";
my $forward_file = "";
my $reverse_file = "";
my $interleaved_file = "";
my $param_file = "";
my $out_dir = "";
my $command;
my $freq = 120;
my $pairend = 0;
my $program = $0;
my $dir = substr($0, 0, length($0)-11);
my $start = time();
my $starttime;
my $endtime;
my $help;
my $assembly_file = "";
my $read_tag = "";
my $mergedFastq = "";
my $max_ext_len = 150;
my $kmin = 0;
my $seq;
my $ncpu ;
my $fgs_complete;
my $fgs_train_reads ;
my $fgs_train_other ;
my $fgs_drop_indel ;
my $assemble ;
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
my $IN1;
my $IN2;
my $val;
my $opt;

#=====================================================================
### Print program usage
##====================================================================
sub print_usage{

  print "\nUSAGE: ./impp_run.pl [options] <fastq file/s> -o <output_directory>";
  print"\n\nInput data options:\n\n";
  print " --12            <filename>    : fastq file with interlaced forward and reverse paired-end reads\n";
  print " -1              <filename>    : fastq file with forward paired-end reads\n";
  print " -2	         <filename>    : fastq file with reverse paired-end reads\n";
  print " -s	         <filename>    : fastq file with unpaired reads\n";
  print " -o/--outdir     <dirname>     : [required] output directory name including the full path\n";
  print " -n              <int>         : [optional] number of threads\n";
  print " -m/--max-len    <int>         : [optional] maximum extension length for anchors [default: 150] \n";
  print " -p/--param-file <filename>    : [optional] parameter file\n";
  print " -h/--help                     : print help message\n";
  print "\n";
  print "NOTE: Input FASTQ read file must be specified in either --12, -1 & -2 or -s format. \n\n"
}


# Getting program options
my($reads_fa, $reads_fq);
getOpts();

#Calling FGS on reads
callFGS_Reads();

#Calling SGA/SPAdes on reads -- adding the mergeSG function to this call
callAssembly();

# Calling mergeSG function -- testing
callMergeSG();

##Split edges file
splitEdges();

##call FGS on split edges
callFGS_SplitEdges();

# Running impp to find and extend anchors
callIMPP();

##Split paths based on length >120bp
splitPaths();

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
print "\nSuccessfully completed running iMPP in : ";
getElapsedTime($endtime - $start);
print "\n";

#====================================================================
## Process the input parameter options
#====================================================================
sub getOpts
{
  my %opts;
  GetOptions(
  \%opts,
  's=s'             =>    \$single_file,
  '1=s'             =>    \$forward_file,
  '2=s'             =>    \$reverse_file,
  '12=s'            =>    \$interleaved_file,
  'outdir|o=s'      =>    \$out_dir,
  'n=s'             =>    \$ncpu,
  'max-len|m=s'     =>    \$max_ext_len,
  'param-file|p=s'  =>    \$param_file,
  'help|h'          =>    \$help
  );

  if ( $help ) {
    print_usage();
    exit;
  }

  if (length($out_dir)==0 || (length($single_file)==0 && length($forward_file)==0 && length($reverse_file)==0 && length($interleaved_file)==0)){
    print "\nERROR: Input parameters were not specified properly. Please see usage below. \n\n";
    print_usage();
    exit;
  }
  else{
    if (length($out_dir) == 0 ){
      print "\nERROR: An output directory name must be specified.\n\n";
      print_usage();
      exit;
    }
    else{
      if(! -e $out_dir)
      {
        $command = "mkdir $out_dir";
        system("$command") == 0 or die "\nError: Could not creat output directory. Please specify full path with -o parameter.\n";
      }
    }
    if (! -e $single_file){
      if(! -e $interleaved_file){
        if( (! -e $forward_file) && (! -e $reverse_file) ){
          print "\nERROR: Input file(s) does not exist. Please try again with a valid input. \n";
          print_usage();
          exit;
        }
        elsif((! -e $forward_file) || (! -e $reverse_file)){
          print "\nERROR: Both forward and reverse reads must be specified. See --help for details.\n\n";
          print_usage();
          exit;
        }
        else{
          if((! -z $forward_file) &&  (! -z $reverse_file)){
            if( ($forward_file=~ /\.fq$/i || $forward_file=~ /\.fastq$/i) && ($reverse_file=~ /(.*).fq/ || $reverse_file=~ /(.*).fastq/)){
              my($IN1, $IN2);
              mergeFastqFiles();
              $pairend = 1;
            }
            else{
              die "\nERROR: Input files must end in .fq or .fastq extension. Please unzip any zipped input files and try again.\n\n";
            }
          }
          else{
            die "\nERROR: Input file(s) are either empty or in compressed format. Please ensure that the reads are in FASTQ format (.fq or .fastq extension) and try again.\n\n";
          }


        }
      }
      else{
        if(! -z $interleaved_file){
          if( ($interleaved_file=~ /\.fq$/i || $interleaved_file=~ /\.fastq$/i)){
            $command = "python ".$dir."utils/fastq2fasta.py $interleaved_file $out_dir";
            system( "$command" ) == 0 or die "\nError: Failed to convert input fastq to fasta file: $? \n";
            $pairend = 1;
          }
          else{
            die "\nERROR: Input files must end in .fq or .fastq extension. Please unzip any zipped input files and try again.\n\n";
          }
        }
        else{
          die "\nERROR: Input file(s) are empty. Please ensure that the input contains reads in FASTQ format and try again.\n\n";
        }


      }
    }
    else{
      if(! -z $single_file){
        if(($single_file=~ /\.fq$/i || $single_file=~ /\.fastq$/i)){
           $command = "python ".$dir."utils/fastq2fasta.py $single_file $out_dir";
           system( "$command" ) == 0 or die "Error: Failed to convert input fastq to fasta file: $?";
        }
        else {  die "\nERROR: Input files must end in .fq or .fastq extension. Please unzip any zipped input files and try again.\n\n";}
      }
      else    { die "\nERROR: Input file(s) are empty. Please ensure that the input contains reads in FASTQ format and try again.\n\n";}
    }


    if (length($assembly_type)==0){ $assembly_type = "0";	}
    #else{
    #    unless ($assembly_type eq "0" || $assembly_type eq "1"){
    #    	print "\nERROR: Please enter a valid input assembler. 0 for sga, 1 for spades. \n\n";
    #    	print_usage();
    #    	exit;
    #	}
    #}

    if ($max_ext_len && $max_ext_len >= 350 ){
      print "\nERROR: Maximum extend length threshold too high. Please select a lower value for optimal results. [default=150] \n\n";
      print_usage();
      exit;
    }

    if (length($param_file) == 0 ){
      if(! -e $dir."params/parameters.txt" ){
        print "\nERROR: Parameter file missing. Make sure parameters.txt file exists in $dir/params/ directory.\n";
        print_usage();
        exit;
      }
      else {  $param_file = $dir."params/parameters.txt"; }
    }
  }

  #Set modified read files path
  $reads_fa = $out_dir."/reads.fasta";
  $reads_fq = $out_dir."/reads.fastq";

  #Loading parameters from parameter file
  loadParameters();

  #Checking if paramters are entered correctly
  checkParameters();
}


#====================================================================
# Call ORFs on input reads using FGS
#====================================================================
sub callFGS_Reads
{
  my $starttime = time();
  print STDERR "\nPredicting ORF's on input reads ...\n";
  print STDERR "\n-----------------------------------------\n";
  $command = $dir."lib/FragGeneScan1.31/run_FragGeneScan.pl";
  $command .= " -genome=$reads_fa";
  $command .= " -out=$out_dir/reads";
  $command .= " -complete=$fgs_complete" ;
  $command .= " -train=$fgs_train_reads";
  $command .= " -thread=$ncpu";
  print "$command\n";
  system( "$command" ) == 0 or die "Error: Could not call FragGeneScan on reads: $?";

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
    if($sga_preprocess){
      $command = "$SGA_BIN preprocess -o assembly.pp.fq $reads_fq";
      $command .= " -p 2" if $pairend;
      system( "$command" ) == 0 or die "Error: Failed to run sga preprocess: $?";
    }
    if($sga_index){
      $command = "$SGA_BIN index -t $ncpu --no-reverse assembly.pp.fq";
      $command .= " -a ropebwt" if $sga_index_algorithm;
      system( "$command" ) == 0 or die "Error: Failed to run sga index: $?";
    }
    if($sga_correct){
      $command = "$SGA_BIN correct -t $ncpu -k $sga_correct_kmer --learn -m $sga_correct_min_olap assembly.pp.fq -o assembly.ec.fq";
      if ($sga_correct_algo == 0) { $command .= " -a kmer";}
      elsif($sga_correct_algo == 1) { $command .= " -a hybrid";}
      elsif($sga_correct_algo == 2) { $command .= " -a overlap";}
      system( "$command" ) == 0 or die "Error: Failed to run sga correct: $?";

      ##index the correct assembly file
      $command = "$SGA_BIN index -t $ncpu -a ropebwt assembly.ec.fq";
      system( "$command" ) == 0 or die "Error: Failed to run sga index on corrected reads: $?";
    }

    if($sga_filter){
      $command = "$SGA_BIN -k $sga_filter_kmer -x $sga_filter_kmer_thresh -t $ncpu assembly.fq ";
      system( "$command" ) == 0 or die "Error: Failed to run sga filter: $?";
    }
    if($sga_rmdup){
      $command = "$SGA_BIN rmdup -t $ncpu assembly.ec.fq";
      system( "$command" ) == 0 or die "Error: Failed to run sga rmdup: $?";
    }
    if($sga_overlap){
      $command = "$SGA_BIN overlap -t $ncpu -m $sga_overlap_min_len assembly.ec.rmdup.fa";
      system( "$command" ) == 0 or die "Error: Failed to run sga overlap: $?";
    }

    #if($sga_assemble){ ## need to remove this for final package
    #$command = "$SGA_BIN assemble -t $ncpu -m $sga_assemble_min_len -g 0 assembly.ec.rmdup.asqg.gz";
    #system( "$command" ) == 0 or die "Error: Failed to run sga assemble: $?";
  #}

  $command = "mv assembly.* $out_dir";
  system($command);

  ## Generating String Graph from the overlap graph output by SGA
  my $olap_file = $out_dir."/assembly.ec.rmdup.asqg.gz";
  print $olap_file;
  if (-e $olap_file) {
    $command = "gunzip -f $olap_file";
    system( "$command" ) == 0 or die "Error: Failed to unzip sga overlap graph file: $?";
  }
  else {
    print "\n Cannot open $olap_file";
  }

  $command = $dir."bin/impp_getSG";
  $command .= " -g $out_dir/assembly.ec.rmdup.asqg -r $out_dir/reads.fasta -f $out_dir/reads.gff";
  system( "$command" ) == 0 or die "Error: Failed to get sga string graph from overlap graph: $?";

  $command = "mv $out_dir/assembly.ec.rmdup.StringGraph.fq  $out_dir/og.fq";
  system($command);
  $assembly_file = $out_dir."/og.fq";

  $endtime = time();

  print STDERR "\n-----------------------------------------\n";
  print STDERR "\nSGA assembly run completed in : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";

  # Calling FGS on assembly edges
  #callFGS_Edges();
  #}
  #elsif ($assembly_type == 1 )
  #{
  ## call SPAdes assembler
  my $SPAdes_bin= $dir."lib/SPAdes-3.9.0-Linux/bin/spades.py";
  $command = "$SPAdes_bin -t $ncpu -m $spades_mem_limit -o $out_dir/spades/assembly";
  if ($pairend) { $command .= " --12 $reads_fq"; }
  else {  $command .= " -s $reads_fq" ;}
  $command .= " -k $spades_kmer";
  $command .= " --only-assembler" if $spades_only_assembler && !$spades_error_correct;
  $command .= " --only-error-correction" if $spades_error_correct && !$spades_only_assembler;
  $command .= " --continue" if $spades_continue;
  $command .= " --disable-gzip-output" if $spades_disable_gzip;
  $command .= " --disable-rr" if $spades_rr;
  system( "$command" ) == 0 or die "Error: Failed to run SPAdes assembler: $?";

  ## modifying the SPAdes graph
  $command = "python ";
  $command .= $dir."utils/renameDB.py";
  $command .= " $out_dir/spades/assembly/assembly_graph.fastg";
  system( "$command" ) == 0 or die "Error: Failed to run modify de bruijn graph file: $?";
  #$assembly_file = $out_dir."/spades/assembly/assembly_graph.dbGraph.fq";

  $endtime = time();
  print STDERR "\n-----------------------------------------\n";
  print STDERR "SPAdes assembly run completed in : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";

  ## Getting kmer Length
  open( DB, '<', $assembly_file ) or die $!;
  while( my $line = <DB> ){
    chomp $line;
    if ( $line =~ /^(>.*)$/ ){
      my $header = $line;
    }
    elsif ( $line !~ /^\s*$/ ){
      $seq = $line;
      chomp $seq;
      my $lseq = length($seq);
      if($kmin == 0)  { $kmin = $lseq;}
      elsif( $lseq < $kmin) { $kmin = $lseq;}
    }
  }

  ## Running bwa on SPAdes contigs
  my $starttime = time();
  print STDERR "\nRunning BWA index (spades contig)...\n";
  $command = $dir."lib/bwa-0.7.17/"."bwa index";
  $command .= " -p $out_dir/spades/contigs.index";
  $command .= " $out_dir/spades/assembly/contigs.fasta";
  print "$command\n";
  system( "$command" ) == 0 or die "Error: Failed to run bwa index on spades contig: $?";


  print STDERR "\nRunning BWA mem...\n";
  $command = $dir."lib/bwa-0.7.17/"."bwa mem";
  $command .= " $out_dir/spades/contigs.index";
  $command .= " $assembly_file";
  $command .= " -p" if $pairend;
  $command .= " -T 45";
  $command .= " -t $ncpu -k $bwa_min_seed"; ## removed -a option
  $command .= " -o $out_dir/og.sam";
  print "$command\n";
  system( "$command" ) == 0 or die "Error: Failed to run bwa mem on edges: $?";
  $endtime = time();
  print "Time taken to complete bwa run on spades contig: \t";
  getElapsedTime($endtime - $starttime);
  print "\n";


}
else{
  print "ERROR: Input assembler value is invalid. Please select either 0 (for sga) or 1 (for spades) .\n"; }

}


#====================================================================
## Call MergeSG function  ( String graph or de Bruijn Graph)
##====================================================================
sub callMergeSG
{
  $assembly_file = $out_dir."/og.fq";


  $starttime = time();
  print STDERR "\nMerging string and de bruijn graph..\n";

  my $command = $dir."bin/impp_mergeSG";
  $command .= " $assembly_file ";
  $command .= " $out_dir/og.sam";
  $command .= " $out_dir/spades/assembly/contigs.fasta ";
  #$command .= " ";
  #$command .= " -k $kmin";

  system( "$command" ) == 0 or die "Error: Failed to run iMPP: $?";


  $endtime = time();

  print STDERR "\n-----------------------------------------\n";
  print "Time taken to merge graphs : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";


}


#====================================================================
# Call split files
#====================================================================
sub splitPaths
{

  $starttime = time();
 #print STDERR "\n-----------------------------------------\n";
  print STDERR "\nSplitting path file to long and short ...";
  print STDERR "\n-----------------------------------------\n";
  $command = $dir."utils/splitFiles.py";
  $command .= " $out_dir/paths.fa";
  $command .=" $freq";
  print "$command\n";
  system( "$command" ) == 0 or die "Error: Failed to run split on paths: $?";
  $endtime = time();
  print STDERR "\n-----------------------------------------\n";
  print "split run on paths completed in : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";

}


#====================================================================
# Call ORFs on assembly graph edges using FGS
#====================================================================
sub splitEdges
{
  $starttime = time();
 #print STDERR "\n-----------------------------------------\n";
  print STDERR "\nSplitting edge file to long and short ...";
  print STDERR "\n-----------------------------------------\n";
  $command = $dir."utils/splitFiles.py";
  $command .= " $out_dir/og.merged.fq";
  $command .=" $freq";
  print "$command\n";
  system( "$command" ) == 0 or die "Error: Failed to run split on assembly edges: $?";
  $endtime = time();
  print STDERR "\n-----------------------------------------\n";
  print "split run on edges completed in : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";

}

#====================================================================
# Call ORFs on assembly graph edges using FGS-split
#====================================================================
sub callFGS_SplitEdges
{
  $starttime = time();
 #print STDERR "\n-----------------------------------------\n";
  print STDERR "\nRunning FragGeneScan on edges ...";
  print STDERR "\n-----------------------------------------\n";
  $command = $dir."lib/FragGeneScan1.31/"."run_FragGeneScan.pl";
  $command .= " -genome=$out_dir/og.merged.short.$freq.fq";
  $command .= " -out=$out_dir/edges.short.$freq";
  $command .= " -complete=0" ;
  $command .= " -train=$fgs_train_other";
  $command .= " -thread=$ncpu";
  print "$command\n";
  system( "$command" ) == 0 or die "Error: Failed to run FragGeneScan on assembly edges: $?";

  $endtime = time();
  print STDERR "\n-----------------------------------------\n";
  print "FragGeneScan run on edges completed in : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";

  $starttime = time();
 #print STDERR "\n-----------------------------------------\n";
  print STDERR "\nRunning FragGeneScan on edges ...";
  print STDERR "\n-----------------------------------------\n";
  $command = $dir."lib/FragGeneScan1.31/"."run_FragGeneScan.pl";
  $command .= " -genome=$out_dir/og.merged.long.$freq.fq";
  $command .= " -out=$out_dir/edges.long.$freq";
  $command .= " -complete=1" ;
  $command .= " -train=complete";
  $command .= " -thread=$ncpu";
  print "$command\n";
  system( "$command" ) == 0 or die "Error: Failed to run FragGeneScan on assembly edges: $?";

  $endtime = time();
  print STDERR "\n-----------------------------------------\n";
  print "FragGeneScan run on edges completed in : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";

  $command = "cat $out_dir/edges.long.$freq.out $out_dir/edges.short.$freq.out > $out_dir/edges.merged.$freq.out";
  system($command);
  $command = "cat $out_dir/edges.long.$freq.gff $out_dir/edges.short.$freq.gff > $out_dir/edges.merged.$freq.gff";
  system($command);
  $command = "cat $out_dir/edges.long.$freq.ffn $out_dir/edges.short.$freq.ffn > $out_dir/edges.merged.$freq.ffn";
  system($command);
  $command = "cat $out_dir/edges.long.$freq.faa $out_dir/edges.short.$freq.faa > $out_dir/edges.merged.$freq.faa";
  system($command);
}



#=====================================================================
## Call ORFs on paths  using FGS
##====================================================================
sub callFGS_Paths
{
  $starttime = time();
  #print STDERR "\n-----------------------------------------\n";
  print STDERR "\nRunning FragGeneScan on candidate paths generated by iMPP ...";
  print STDERR "\n-----------------------------------------\n";
  $command = $dir."lib/FragGeneScan1.31/"."run_FragGeneScan.pl";
  $command .= " -genome=$out_dir/paths.short.$freq.fq";
  $command .= " -out=$out_dir/paths.short.$freq";
  $command .= " -complete=0" ;
  $command .= " -train=$fgs_train_other";
  $command .= " -thread=$ncpu";
  print "$command\n";
  system( "$command" ) == 0 or die "Error: Failed to run FragGeneScan on candidate paths: $?";

  $endtime = time();
  print STDERR "\n-----------------------------------------\n";
  print "FragGeneScan run on paths completed in : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";

  $starttime = time();
  #print STDERR "\n-----------------------------------------\n";
  print STDERR "\nRunning FragGeneScan on candidate paths generated by iMPP ...";
  print STDERR "\n-----------------------------------------\n";
  $command = $dir."lib/FragGeneScan1.31/"."run_FragGeneScan.pl";
  $command .= " -genome=$out_dir/paths.long.$freq.fq";
  $command .= " -out=$out_dir/paths.long.$freq";
  $command .= " -complete=1" ;
  $command .= " -train=complete";
  $command .= " -thread=$ncpu";
  print "$command\n";
  system( "$command" ) == 0 or die "Error: Failed to run FragGeneScan on candidate paths: $?";

  $endtime = time();
  print STDERR "\n-----------------------------------------\n";
  print "FragGeneScan run on paths completed in : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";

  $command = "cat $out_dir/paths.long.$freq.out $out_dir/paths.short.$freq.out > $out_dir/paths.merged.$freq.out";
  system($command);
  $command = "cat $out_dir/paths.long.$freq.gff $out_dir/paths.short.$freq.gff > $out_dir/paths.merged.$freq.gff";
  system($command);
  $command = "cat $out_dir/paths.long.$freq.ffn $out_dir/paths.short.$freq.ffn > $out_dir/paths.merged.$freq.ffn";
  system($command);
  $command = "cat $out_dir/paths.long.$freq.faa $out_dir/paths.short.$freq.faa > $out_dir/paths.merged.$freq.faa";
  system($command);
}


#======================================================================
### Aligning reads to edges and paths using BWA
###====================================================================
sub callBWA
{
  return unless $bwa;
  my $starttime = time();
  print STDERR "\nRunning BWA index (edges)...\n";
  $command = $dir."lib/bwa-0.7.17/"."bwa index";
  $command .= " -p $out_dir/edges.index";
  $command .= " $out_dir/og.merged.fq";
  print "$command\n";
  system( "$command" ) == 0 or die "Error: Failed to run bwa index on edges: $?";

  print STDERR "\nRunning BWA mem...\n";
  $command = $dir."lib/bwa-0.7.17/"."bwa mem";
  $command .= " $out_dir/edges.index";
  $command .= " $reads_fq";
  $command .= " -p" if $pairend;
  $command .= " -T $bwa_min_score";
  $command .= " -a -t $ncpu -k $bwa_min_seed"; # removed -a option
  $command .= " -o $out_dir/edges.$bwa_min_score.sam";
  print "$command\n";
  system( "$command" ) == 0 or die "Error: Failed to run bwa mem on edges: $?";
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
  system( "$command" ) == 0 or die "Error: Failed to run bwa index on candidate paths: $?";

  print STDERR "\nRunning BWA mem...\n";
  $command = $dir."lib/bwa-0.7.17/"."bwa mem";
  $command .= " $out_dir/paths.index";
  $command .= " $reads_fq";
  $command .= " -p" if $pairend;
  $command .= " -T $bwa_min_score";
  $command .= " -a -t $ncpu -k $bwa_min_seed"; ##removed -a option
  $command .= " -o $out_dir/paths.$bwa_min_score.sam";
  print "$command\n";
  system( "$command" ) == 0 or die "Error: Failed to run bwa mem on candidate paths: $?";

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
  print STDERR "\nExtending anchors..\n";

  my $command = $dir."bin/impp_extendAnchor";
  $command .= " -g $out_dir/og.merged.fq";
  $command .= " -a $out_dir/edges.merged.$freq.gff";
  $command .= " -l $max_ext_len ";
  $command .= " -p $out_dir/paths.fa";
  $command .= " -k $kmin";

  system( "$command" ) == 0 or die "Error: Failed to run iMPP: $?";


  $endtime = time();

  print STDERR "\n-----------------------------------------\n";
  print "Time taken to find and extend Anchors : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";
}


#=====================================================================
## Generating six frame translation
###====================================================================
sub callSixFrameFilter
{
  $starttime = time();
  #print STDERR "\n-----------------------------------------\n";
  print STDERR "\nRunning six frame translation filter on any un-called reads..";
  print STDERR "\n---------------------------------------------------------------\n";
  my $command = "python3 ";
  $command .= $dir."utils/sixFrameTranslate.py";
  $command .= " $out_dir/reads.gff ";
  $command .= " $out_dir/edges.merged.$freq.gff ";
  $command .= " $out_dir/paths.merged.$freq.gff ";
  $command .= " $out_dir/edges.$bwa_min_score.sam";
  $command .= " $out_dir/paths.$bwa_min_score.sam ";
  $command .= " $reads_fa";
  $command .= " $out_dir";
  system( "$command" ) == 0 or die "Error: Failed to run iMPP six frame filter: $?";
  $endtime = time();

  #print STDERR "\n-----------------------------------------\n";
  print "Time taken to generate six frame translation : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";
  $command = "mv $out_dir/reads.filter.fasta $out_dir/reads.filter.$bwa_min_score.fasta";
  system($command);
  $command = "mv $out_dir/reads.translated.fasta $out_dir/reads.translated.$bwa_min_score.fasta";
  system($command);
}

#=====================================================================
### Generating final peptide assembly
####====================================================================
sub callPeptideAssembler
{

  return unless $plass;
  $starttime = time();
  print STDERR "\nGenerating peptide assembly..\n";
  print STDERR "\n-----------------------------------------\n";
  if ( -e $out_dir."/assembled_proteins.$bwa_min_score.faa")
  {
    $command = "rm -f $out_dir/assembled_proteins.$bwa_min_score.faa";
    system($command);
  }
  if( -d $out_dir."/plass.tmp")
  {
    $command = "rm -rf $out_dir/plass.tmp";
    system($command);
  }
  $command = $dir."lib/plass/bin/plass";
  $command .= " assemble" if $plass_assemble;
  $command .= " --threads $ncpu --min-length $plass_min_len --num-iterations $plass_num_iter";
  $command .= " $out_dir/reads.filter.$bwa_min_score.fasta $out_dir/assembled_proteins.$bwa_min_score.faa";
  $command .= " $out_dir/plass.tmp ";

  system( "$command" ) == 0 or die "Error: Failed to run peptide assembly: $?";
  system("rm -rf $out_dir/plass.tmp");

  $endtime = time();

  print STDERR "\n-----------------------------------------\n";
  print "Time taken to complete peptide assembly : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";

  print STDERR "\nFiltering plass generated contigs..\n";
  print STDERR "\n-----------------------------------------\n";

  $command = "python ";
  $command .= $dir."utils/renamePlassContig.py";
  $command .= " $out_dir/assembled_proteins.$bwa_min_score.faa";
  system( "$command" ) == 0 or die "Error: Failed to filter plass contigs: $?";

  ##getting mapping of six frame translated reads to plass contigs -- using diamond Blastx
  $command = $dir."lib/diamond";
  $command .= " makedb --in $out_dir/assembled_proteins.$bwa_min_score.re.60.faa -d $out_dir/plass_contig.index";
  system( "$command" ) == 0 or die "Error: Failed to generate index for plass conitigs: $?";
  my $command = $dir."lib/diamond";
  $command .= " blastx -p $ncpu -d $out_dir/plass_contig.index -q $out_dir/reads.translated.$bwa_min_score.fasta -o $out_dir/plass.translated.$bwa_min_score.dmd";
  system( "$command" ) == 0 or die "Error: Failed to run blastx on reads against plass contigs: $?";



  ##getting mapped translated reads from diamond blastx mapping
  $command = "python ";
  $command .= $dir."utils/plassMapping.py";
  $command .= " $out_dir/assembled_proteins.$bwa_min_score.re.60.faa";
  $command .= " $out_dir/plass.translated.$bwa_min_score.dmd";
  $command .= " $out_dir/reads.translated.$bwa_min_score.fasta";
  system( "$command" ) == 0 or die "Error: Failed to filter plass contigs (mapped reads): $?";

  ##getting impp refined reads for second round of plass run
  $command = "python ";
  $command .= $dir."utils/getPredReads.py";
  $command .= " $out_dir/reads.filter.$bwa_min_score.fasta";
  $command .= " $out_dir/reads.translated.$bwa_min_score.fasta";
  system( "$command" ) == 0 or die "Error: Failed to filter plass contigs (mapped reads): $?";
  $command = "cat $out_dir/reads.filter.$bwa_min_score.pred.fasta $out_dir/reads.translated.$bwa_min_score.mapped.fasta > $out_dir/reads.impp.fasta";
  system( "$command" ) == 0 or die "Error: Failed to get impp filtered reads: $?";

  ##run plass with impp reads 
  $command = $dir."lib/plass/bin/plass";
  $command .= " assemble" if $plass_assemble;
  $command .= " --threads $ncpu --min-length $plass_min_len --num-iterations $plass_num_iter";
  $command .= " $out_dir/reads.impp.fasta $out_dir/assembled_proteins.$bwa_min_score.impp.faa";
  $command .= " $out_dir/plass.tmp ";
  system( "$command" ) == 0 or die "Error: Failed to run peptide assembly: $?";
  system("rm -rf $out_dir/plass.tmp");
  
  $endtime = time();
  print STDERR "\n-----------------------------------------\n";
  print "Time taken to complete peptide contigs filtering : ";
  getElapsedTime($endtime - $starttime);
  print STDERR "\n-----------------------------------------\n";

}


#==================================================
# Program to merge seperate forward and reverse paired end read files into single interleaved file
sub mergeFastqFiles
{
  open ($IN1,'<', $forward_file);
  open ($IN2, '<' ,$reverse_file);
  open (my $mergedFastq, '>', $out_dir."/reads.merged.fq");
  while(<$IN1>) {
    print $mergedFastq $_;
    $_ = <$IN1>;
    print $mergedFastq $_;
    $_ = <$IN1>;
    print $mergedFastq $_;
    $_ = <$IN1>;
    print $mergedFastq $_;

    $_ = <$IN2>;
    print $mergedFastq $_;
    $_ = <$IN2>;
    print $mergedFastq $_;
    $_ = <$IN2>;
    print $mergedFastq $_;
    $_ = <$IN2>;
    print $mergedFastq $_;
  }

  close($IN1);
  close($IN2);

  #Program to convert fastq to fasta
  $command = "python ".$dir."utils/fastq2fasta.py ".$out_dir."/reads.merged.fq $out_dir";
  print $command;
  system( "$command" ) == 0 or die "Error: Failed to convert input fastq to fasta: $?";
  $command = "rm -f".$out_dir."/reads.merged.fq";

}

#=====================================================================
#### Loading parameters for running 3rd party softwares
#####=================================================================
sub loadParameters
{
  open(IN, '<', $param_file) or die $!;
  while ( <IN>) {
    chomp;
    next if /^#/;
    next unless /\S/;
    ($opt, $val)  = split();
    for ( $opt ) {
      if ($opt eq "thread")		                 {$ncpu     = $val; }
      if ($opt eq "mlen")                      {$max_ext_len = $val;}
      if ($opt eq "fgs_complete")              {$fgs_complete = $val; }
      if ($opt eq "fgs_train_reads")           {$fgs_train_reads    = $val; }
      if ($opt eq "fgs_train_other")           {$fgs_train_other    = $val; }
      if ($opt eq "fgs_dropindel")             {$fgs_drop_indel   = $val; }
      if ($opt eq "assemble")                  {$assemble   = $val; }
      if ($opt eq "sga_preprocess")             {$sga_preprocess = $val;}
      if ($opt eq "sga_index")                  {$sga_index = $val;}
      if ($opt eq "sga_index_algorithm")       {$sga_index_algorithm = $val;}
      if ($opt eq "sga_overlap")                {$sga_overlap = $val;}
      if ($opt eq "sga_overlap_min_len")        {$sga_overlap_min_len = $val;}
      if ($opt eq "sga_rmdup")                  {$sga_rmdup = $val;}
      if ($opt eq "sga_assemble")               {$sga_assemble = $val;}
      if ($opt eq "sga_assemble_min_len")       {$sga_assemble_min_len = $val;}
      if ($opt eq "sga_merge")                  {$sga_merge = $val;}
      if ($opt eq "sga_bwt2fa")                 {$sga_bwt2fa = $val;}
      if ($opt eq "sga_correct")                {$sga_correct = $val;}
      if ($opt eq "sga_correct_algo")           {$sga_correct_algo = $val;}
      if ($opt eq "sga_correct_kmer")           {$sga_correct_kmer = $val;}
      if ($opt eq "sga_correct_min_olap")       {$sga_correct_min_olap = $val;}
      if ($opt eq "sga_fmmerge")                {$sga_fmmerge = $val;}
      if ($opt eq "sga_fmmerge_min_olap")       {$sga_fmmerge_min_olap = $val;}
      if ($opt eq "sga_filter")                 {$sga_filter = $val;}
      if ($opt eq "sga_filter_kmer")            {$sga_filter_kmer = $val;}
      if ($opt eq "sga_filter_kmer_thresh")     {$sga_filter_kmer_thresh = $val;}
      if ($opt eq "spades")                     {$spades = $val;}
      if ($opt eq "spades_only_assembler")      {$spades_only_assembler = $val;}
      if ($opt eq "spades_error_correct")       {$spades_error_correct = $val;}
      if ($opt eq "spades_memory_limit")        {$spades_mem_limit = $val;}
      if ($opt eq "spades_disable_gzip")        {$spades_disable_gzip = $val;}
      if ($opt eq "spades_continue")            {$spades_continue = $val;}
      if ($opt eq "spades_kmer")                {$spades_kmer = $val;}
      if ($opt eq "spades_rr")                  {$spades_rr = $val;}
      if ($opt eq "bwa")                        {$bwa = $val;}
      if ($opt eq "bwa_min_seed")               {$bwa_min_seed = $val;}
      if ($opt eq "bwa_min_score")              {$bwa_min_score = $val;}
      if ($opt eq "plass")                      {$plass = $val;}
      if ($opt eq "plass_assemble")             {$plass_assemble = $val;}
      if ($opt eq "plass_num_iter")             {$plass_num_iter = $val;}
      if ($opt eq "plass_min_length")           {$plass_min_len = $val;}

    }
  }
  close(IN);

}

#=====================================================================
#### Program to check input parameters in parameter file
#####=================================================================

sub checkParameters
{
  my $max_cpu = `grep -w "processor" /proc/cpuinfo | wc -l`; chomp $max_cpu;
  if ( $ncpu > $max_cpu ) {$ncpu = $max_cpu; }
  if ( ! defined $max_ext_len ) { die "ERROR: Maximum extend length must be given\n"; }
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
  if ( $bwa_min_score < 0)  { die "ERROR: Invalid bwa min score option\n"; }
  if ( $plass != 0 && $plass != 1 ) { die "ERROR: Invalid plass option\n"; }
  if ( $plass_assemble != 0 && $plass_assemble != 1 ) { die "ERROR: Invalid plass assemble option\n"; }
  if ( $plass_min_len < 0 ) { die "ERROR: Invalid plass min length option\n"; }
  if ( $plass_num_iter < 0 ) { die "ERROR: Invalid plass num iterations option\n"; }


}


#=====================================================================
### Clean up final output directory
###====================================================================
sub cleanOutput
{
  my $plass_dir = "$out_dir/plass.tmp/";
  my $spades_dir = "$out_dir/spades/";
  $command = "cat $out_dir/edges.merged.$freq.gff $out_dir/paths.merged.$freq.gff > $out_dir/orfs.gff";
  system($command);
  $command = "cat $out_dir/edges.merged.$freq.ffn $out_dir/paths.merged.$freq.ffn > $out_dir/orfs.ffn";
  system($command);
  $command = "cat $out_dir/edges.merged.$freq.faa $out_dir/paths.merged.$freq.faa > $out_dir/orfs.faa";
  system($command);
  $command = "cat $out_dir/assembled_proteins.$bwa_min_score.faa > $out_dir/assembled_proteins.faa";
  system($command);
  unlink glob "'$out_dir/og.*'";
  unlink glob "'$out_dir/reads*.*'";
  unlink glob "'$out_dir/edges*.*'";
  unlink glob "'$out_dir/assembly*.*'";
  unlink glob "'$out_dir/paths*.*'";
  unlink glob "'$out_dir/*dmnd'";
  unlink glob "'$out_dir/*dmd'";
  unlink glob "'$out_dir/assembled_proteins.$bwa_min_score.*'";

  rmtree $spades_dir;
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
