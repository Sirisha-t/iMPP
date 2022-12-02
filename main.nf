#! /usr/bin/env nextflow

nextflow.enable.dsl=2


/** Print help message **/
def helpMessage(){
log.info """
 USAGE: nextflow run main.nf --single/--forward --reverse/--interleaved [other options] <fastq file/s> --outdir <output_directory> -profile docker
 Input data options:
   -profile               <string>      : [required] docker, base/test
   --interleaved          <filename>    : fastq file with interlaced forward and reverse paired-end reads
   --forward              <filename>    : fastq file with forward paired-end reads
   --reverse              <filename>    : fastq file with reverse paired-end reads
   --single               <filename>    : fastq file with unpaired reads
   --outdir               <dirname>     : [required] output directory name including the full path [default: output]
   --genecaller           <string>      : [optional] both, fgs or prodigal [default: both]
   --threads              <int>         : [optional] number of threads [default: 16]
   --maxlen               <int>         : [optional] maximum extension length for anchors [default: 150]
   -h/--help                            :  help message

  NOTE: Input FASTQ read file must be specified in either --interleaved, --forward & --reverse or --single format.

"""
}

fgs = false
mpd = false

// Show help message
if (params.help){
  helpMessage()
  exit 0
}

/** Checking if input parameters have been entered correctly **/
if (params.single == "null" && params.forward=="null" && params.reverse=="null" && params.interleaved=="null"){
    println("\nERROR: Input fastq files were not specified properly. Please see usage below. \n\n")
    helpMessage()
    exit 1
}
else{
if(params.single != "null") paired = 0
if(params.reverse != "null" && params.forward != "null")  paired = 1
if(params.interleaved != "null")  paired =1
}


if(params.outdir == "null" ){
    println("\nERROR: OUtput dir not specified. Please see usage below. \n\n")
    helpMessage()
    exit 1
}
else{
    outpath = file(params.outdir)
  if(!outpath.exists())
  {
    result = outpath.mkdir()
    println result ? "Created output dir" : "Cannot create directory: $outpath"
  }

}


if (params.maxlen && params.maxlen >= 350 ){
      println("\nERROR: Maximum extend length threshold too high. Please select a lower value for optimal results. [default=150] \n\n")
      helpMessage()
      exit 1
    }

if(params.genecaller == "both"){
    fgs = true
    mpd = true
  }
else if(params.genecaller == "fgs") fgs = true
else if(params.genecaller == "prodigal")  mpd = true
else{
    println("ERROR: Invalid gene caller option specified. Please use either both, fgs or prodigal. [default=both]")
    helpMessage()
    exit 1
  }


workflow FGS_EDGE{
    take:
    reads_gff
    sga_assembly
    spades_assembly
    reads_fastq
    reads_fasta

    main:
    getStringGraph(sga_assembly, reads_fastq, reads_gff, params.fgs_tag)
    callBWAContigs(spades_assembly, getStringGraph.out.string_graph, params.fgs_tag)
    getMergedGraph(getStringGraph.out.string_graph,callBWAContigs.out.bwa_spades_contig, spades_assembly, params.fgs_tag)
    splitInput(getMergedGraph.out.merged_graph, params.fgs_tag)
    callFGS(splitInput.out.short_file, splitInput.out.long_file, params.edge_out, params.fgs_tag)

    FGS_PATH(getMergedGraph.out.merged_graph, callFGS.out.fgs_gff, reads_fastq, reads_fasta, reads_gff, )
    
}


workflow FGS_PATH{
    take:
    merged_graph
    fgs_edges
    reads_fastq
    reads_fasta
    reads_gff
    
    main:
    callIMPP(merged_graph, fgs_edges, params.path_out, params.fgs_tag)
    splitInput(callIMPP.out.impp_paths, params.fgs_tag)
    callFGS(splitInput.out.short_file, splitInput.out.long_file, params.path_out, params.fgs_tag)

    callBWA(merged_graph, callIMPP.out.impp_paths, reads_fastq, paired, params.fgs_tag)
    sixFrameTranslate(reads_fasta, reads_gff, fgs_edges, callFGS.out.fgs_gff, callBWA.out.bwa_edge, callBWA.out.bwa_path, params.fgs_tag)
    peptideAssembly_first(sixFrameTranslate.out.filtered_reads, params.fgs_tag)
    modifyAssemblyOut(peptideAssembly_first.out.plass_first_round, params.fgs_tag)
    alignSixFrametoAssembly(modifyAssemblyOut.out.plass_renamed, sixFrameTranslate.out.translated_reads, params.fgs_tag)
    refinePlassInput(modifyAssemblyOut.out.plass_renamed, alignSixFrametoAssembly.out.plass_alignment, sixFrameTranslate.out.translated_reads, sixFrameTranslate.out.filtered_reads, params.fgs_tag )
    peptideAssembly_second(refinePlassInput.out.impp_filtered_reads, params.fgs_tag)
   
}   
    
workflow MPD_EDGE{
    take:
    reads_gff
    sga_assembly
    spades_assembly
    reads_fastq
    reads_fasta

    main:
    getStringGraph(sga_assembly, reads_fastq, reads_gff, params.mpd_tag)
    callBWAContigs(spades_assembly, getStringGraph.out.string_graph, params.mpd_tag)
    getMergedGraph(getStringGraph.out.string_graph,callBWAContigs.out.bwa_spades_contig, spades_assembly, params.mpd_tag)
    splitInput(getMergedGraph.out.merged_graph, params.mpd_tag)
    callMPD(splitInput.out.short_file, splitInput.out.long_file, params.edge_out, params.mpd_tag)

    MPD_PATH(getMergedGraph.out.merged_graph, callMPD.out.mpd_gff,reads_fastq, reads_fasta, reads_gff)

}

workflow MPD_PATH{
    take:
    merged_graph
    mpd_edges
    reads_fastq
    reads_fasta
    reads_gff
    
    main:
    callIMPP(merged_graph, mpd_edges, params.path_out, params.mpd_tag)
    splitInput(callIMPP.out.impp_paths, params.mpd_tag)
    callMPD(splitInput.out.short_file, splitInput.out.long_file, params.path_out, params.mpd_tag)
    
    callBWA(merged_graph, callIMPP.out.impp_paths, reads_fastq, paired, params.mpd_tag)
    sixFrameTranslate(reads_fasta, reads_gff, mpd_edges, callMPD.out.mpd_gff, callBWA.out.bwa_edge, callBWA.out.bwa_path, params.mpd_tag)
    peptideAssembly_first(sixFrameTranslate.out.filtered_reads, params.mpd_tag)
    modifyAssemblyOut(peptideAssembly_first.out.plass_first_round, params.mpd_tag)
    alignSixFrametoAssembly(modifyAssemblyOut.out.plass_renamed, sixFrameTranslate.out.translated_reads, params.mpd_tag)
    refinePlassInput(modifyAssemblyOut.out.plass_renamed, alignSixFrametoAssembly.out.plass_alignment, sixFrameTranslate.out.translated_reads, sixFrameTranslate.out.filtered_reads, params.mpd_tag )
    peptideAssembly_second(refinePlassInput.out.impp_filtered_reads, params.mpd_tag)
}   


workflow{
    if (params.single == "null"){ 
      if(params.interleaved == "null"){
        if( (params.forward == "null") && (params.reverse == "null") ){
          println("\nERROR: Input file(s) does not exist. Please try again with a valid input. \n")
          helpMessage()
          exit 1
        }
        else if((params.forward == "null") || (params.reverse == "null")){
          println("\nERROR: Both forward and reverse reads must be specified. See --help for details.\n\n")
          helpMessage()
          exit 1
        }
        else{
	        forwardfile = file(params.forward)
	        reversefile = file(params.reverse)
	        if( (forwardfile.getExtension() ==  "fq" || forwardfile.getExtension() ==  "fastq"|| forwardfile.getExtension() == "gz") && (reversefile.getExtension() ==  "fq" || reversefile.getExtension() ==  "fastq" || reversefile.getExtension() == "gz")){
              fwd_file = Channel.fromPath(params.forward)
              rev_file = Channel.fromPath(params.reverse)
              outputfile = Channel.fromPath(outpath)
              mergeFastqFiles(fwd_file, rev_file)
              convertFastq(mergeFastqFiles.out.merged_file, outputfile) 
              params.paired = 1
            }
            else{
              println("\nERROR: Input files must end in .fq or .fastq extension. Please unzip any zipped input files and try again.\n\n")
              exit 1
            }
          }

    }
    else{
	      fastq_file = file(params.interleaved)
          if((fastq_file.getExtension() ==  "fq"|| fastq_file.getExtension() == "fastq" || fastq_file.getExtension() == "gz")){
            inputfile = Channel.fromPath(params.interleaved)
            outputfile = Channel.fromPath(outpath)
            convertFastq(inputfile, outputfile) 
            params.paired = 1
          }
          else{
            println("\nERROR: Input files must end in .fq or .fastq extension. Please unzip any zipped input files and try again.\n\n")
            exit 1
          }
        }
    }


else{
    fastq_file = file(params.single)
      if((fastq_file.getExtension() ==  "fq"|| fastq_file.getExtension() == "fastq" || fastq_file.getExtension() == "gz")){   
          inputfile = Channel.fromPath(params.single)
          outputfile = Channel.fromPath(outpath)
          convertFastq(inputfile, outputfile)  
          params.paired = 0
      }
      else { 
          println("\nERROR: Input files must end in .fq or .fastq extension. Please unzip any zipped input files and try again.\n\n")
          exit 1
      }
}
  
//generate overlap graph using SGA
callSGAAssembly(convertFastq.out.reads_fastq, paired)
modifyASQG(callSGAAssembly.out.asqg_file)

//Generate de Bruijn graph using SPAdes
callSPAdesAssembly(convertFastq.out.reads_fastq, paired)

if( fgs == true && mpd == true)
{
  //run FGS-based pipeline
  callFGSRead(convertFastq.out.reads_fasta, params.fgs_tag)
  FGS_EDGE(callFGSRead.out.reads_gff, modifyASQG.out.asqg_modfile, callSPAdesAssembly.out.contig, convertFastq.out.reads_fastq, convertFastq.out.reads_fasta)
        
  //run MPD-based pipeline  
  callMPDRead(convertFastq.out.reads_fasta, params.mpd_tag) 
  MPD_EDGE(callMPDRead.out.reads_gff,modifyASQG.out.asqg_modfile, callSPAdesAssembly.out.contig, convertFastq.out.reads_fastq, convertFastq.out.reads_fasta)
}
else{
  if(fgs == true){
    callFGSRead(convertFastq.out.reads_fasta, params.fgs_tag)
    FGS_EDGE(callFGSRead.out.reads_gff, modifyASQG.out.asqg_modfile, callSPAdesAssembly.out.contig, convertFastq.out.reads_fastq, convertFastq.out.reads_fasta)
  }
  else{
    callMPDRead(convertFastq.out.reads_fasta, params.mpd_tag) 
    MPD_EDGE(callMPDRead.out.reads_gff,modifyASQG.out.asqg_modfile, callSPAdesAssembly.out.contig, convertFastq.out.reads_fastq, convertFastq.out.reads_fasta)
  }
}


}


process convertFastq{
  cache 'deep'
  //publishDir "${params.tmpdir}"
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path queryFile
  path out

  output:
  path 'reads.fasta', emit: reads_fasta
  path 'reads.fastq', emit: reads_fastq
  """
  sed -n '1~4s/^@/>/p;2~4p' $queryFile > reads.fasta
  cat $queryFile > reads.fastq
  """
}


process mergeFastqFiles{
  cache 'deep'
  //publishDir "${params.tmpdir}"
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path fwd
  path rev

  output:
  path 'merged.fastq', emit: merged_file

  script:
  """
  seqtk mergepe $fwd $rev > merged.fastq
  """
}

/**Call processes based on workflow execution**/

process callFGSRead{
  cache 'deep'
  executor 'local'
  cpus params.threads
  
  input:
    path(queryFile)
    val(tag)
  output:
    //publishDir "${params.tmpdir}"
    publishDir "${params.outdir}", mode: 'copy'
    path "${tag}.reads.gff", emit: reads_gff
    path "${tag}.reads.ffn", emit: reads_ffn
    path "${tag}.reads.faa", emit: reads_faa
    path "${tag}.reads.out", emit: reads_out

  script:
  """
  run_FragGeneScan.pl -thread=$params.threads -genome=$queryFile -out=${tag}.reads -complete=$params.fgs_mode -train=$params.fgs_train
  """
}

process callMPDRead{
  executor 'local'
  cpus params.threads
  
  input:
    path(queryFile)
    val(tag)
  output:
    //publishDir "${params.tmpdir}"
    publishDir "${params.outdir}", mode: 'copy'
    path "${tag}.reads.gff", emit: reads_gff
    path "${tag}.reads.ffn", emit: reads_ffn
    path "${tag}.reads.faa", emit: reads_faa
    path "${tag}.reads.genes", emit: reads_out
  script:
  """
  prodigal -i $queryFile -d ${tag}.reads.ffn -f gff -o ${tag}.reads.gff -p meta -q -s ${tag}.reads.genes -a ${tag}.reads.faa
  """
}


process callSGAAssembly{
  cache 'deep'
  executor 'local'
  cpus params.threads
  

  input:
  path(queryFile)
  val(paired)
  output:
  path 'assembly.*', emit: sga_assembly
  path 'assembly.ec.rmdup.asqg', emit: asqg_file
  //publishDir "${params.tmpdir}"
  publishDir "${params.outdir}", mode: 'copy'
  script:
    if(paired == 1){
    """
    sga preprocess -o assembly.pp.fq -p 2 $queryFile
    sga index --no-reverse -t $params.threads -d $params.sga_disk assembly.pp.fq
    sga correct -t $params.threads -k $params.sga_correct_kmer --learn -m $params.sga_correct_min_olap assembly.pp.fq -o assembly.ec.fq
    sga index -t $params.threads -a ropebwt assembly.ec.fq
    sga rmdup -t $params.threads assembly.ec.fq
    sga overlap -t $params.threads -m $params.sga_overlap_min_len assembly.ec.rmdup.fa
    gunzip -f assembly.ec.rmdup.asqg.gz
    """
    }

    else{
    """
    sga preprocess -o assembly.pp.fq $queryFile
    sga index --no-reverse -t $params.threads -d $params.sga_disk assembly.pp.fq
    sga correct -t $params.threads -k $params.sga_correct_kmer --learn -m $params.sga_correct_min_olap assembly.pp.fq -o assembly.ec.fq
    sga index -t $params.threads -a ropebwt assembly.ec.fq
    sga rmdup -t $params.threads assembly.ec.fq
    sga overlap -t $params.threads -m $params.sga_overlap_min_len assembly.ec.rmdup.fa
    gunzip -f assembly.ec.rmdup.asqg.gz
    """
    }
    

}

process callSPAdesAssembly{
  cache 'deep'
  executor 'local'
  cpus params.threads
  //publishDir "${params.tmpdir}"

  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(queryFile)
  val(paired)
  output:
  path 'spades/*' , emit: spades_outdir
  path 'spades/assembly_graph.dbGraph.*', emit: spades_assembly
  path 'spades/contigs.fasta', emit: contig

  script:
    if(paired == 1)
    {
    """
    spades.py -t $params.threads -m $params.spades_memory_limit -o $params.spadesdir --12 $queryFile --meta --only-assembler -k $params.spades_kmer 
    renameDB.py $params.spadesdir/assembly_graph.fastg
    """
    }
    else{
    """
    spades.py -t $params.threads -m $params.spades_memory_limit -o $params.spadesdir -s $queryFile --only-assembler -k $params.spades_kmer 
    renameDB.py $params.spadesdir/assembly_graph.fastg
    """
    }
  

}

process modifyASQG{
  cache 'deep'
  publishDir "${params.outdir}", mode: 'copy'
  
  input:
  path(asqg)
  
  output:
  path 'assembly.ec.rmdup.mod.asqg', emit: asqg_modfile
  path 'assembly.ec.rmdup.namemap.txt', emit: asqg_namemap
  """
  impp_modASQG -g $asqg
  """
}

process getStringGraph{
  cache 'deep'
  publishDir "${params.outdir}", mode: 'copy'
  
  input:
  path(asqg)
  path(queryFile)
  path(reads_gff)
  val(tag)
  
  output:
  
  path "${tag}.og.fq", emit: string_graph

  script:
    """
    impp_getSG -g $asqg -r $queryFile -f $reads_gff
    mv assembly.ec.rmdup.mod.StringGraph.fq  ${tag}.og.fq
    """
}

process callBWAContigs{
  cache 'deep'
  cpus params.threads 
  container 'biocontainers/bwa:v0.7.17_cv1'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(spades_contig)
  path(string_graph)
  val(tag)
  output:
  path "${tag}.og.sam", emit: bwa_spades_contig
  """
  bwa index -p contigs.index $spades_contig
  bwa mem contigs.index $string_graph -T 45 -t $params.threads -k $params.bwa_min_seed -o ${tag}.og.sam
  """
}


process getMergedGraph{
  cache 'deep'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(string_graph)
  path(bwa_contig)
  path(spades_contig)
  val(tag)
  output:
  path "${tag}.og.merged.fq", emit: merged_graph
  """
  impp_mergeSG $string_graph $bwa_contig $spades_contig
  """
} 

process splitInput{
  cache 'deep'
  container 'staphb/spades:latest'
  //publishDir "${params.tmpdir}"
  publishDir "${params.outdir}", mode: 'copy'
  input:
  path(mergedFile)
  val(tag)
  output:
  path "${mergedFile.baseName}.short.*", emit: short_file
  path "${mergedFile.baseName}.long.*", emit: long_file

  """
  splitFiles.py $mergedFile $params.freq 
  """

}


process callFGS{
  executor 'local'
  cpus params.threads
  container 'tsirisha/fraggenescan:latest'
  //publishDir "${params.tmpdir}"
  publishDir "${params.outdir}", mode: 'copy'
  input:
  path(shortFile)
  path(longFile)
  val(outname)
  val(tag)
  output:
    path "${tag}.${outname}.merged.${params.freq}.gff", emit: fgs_gff
    path "${tag}.${outname}.merged.${params.freq}.ffn", emit: fgs_ffn
    path "${tag}.${outname}.merged.${params.freq}.faa", emit: fgs_faa
    path "${tag}.${outname}.merged.${params.freq}.out", emit: fgs_out
  script:
  """
  
  run_FragGeneScan.pl -genome=${shortFile.baseName}.fq -out=${tag}.${outname}.short.${params.freq} -complete=0 -train=$params.fgs_train_other -thread=$params.threads
  run_FragGeneScan.pl -genome=${longFile.baseName}.fq -out=${tag}.${outname}.long.${params.freq} -complete=1 -train=$params.fgs_train_complete -thread=$params.threads
  cat ${tag}.${outname}.long.${params.freq}.out ${tag}.${outname}.short.${params.freq}.out > ${tag}.${outname}.merged.${params.freq}.out
  cat ${tag}.${outname}.long.${params.freq}.gff ${tag}.${outname}.short.${params.freq}.gff > ${tag}.${outname}.merged.${params.freq}.gff
  cat ${tag}.${outname}.long.${params.freq}.ffn ${tag}.${outname}.short.${params.freq}.ffn > ${tag}.${outname}.merged.${params.freq}.ffn
  cat ${tag}.${outname}.long.${params.freq}.faa ${tag}.${outname}.short.${params.freq}.faa > ${tag}.${outname}.merged.${params.freq}.faa
  """
}

process callMPD{
  executor 'local'
  cpus params.threads
  container 'biocontainers/prodigal:v1-2.6.3-4-deb_cv1'
  publishDir "${params.outdir}", mode: 'copy'
  
  input:
    path(shortFile)
    path(longFile)
    val(outname)
    val(tag)
  output:
    path "${tag}.${outname}.merged.${params.freq}.gff", emit: mpd_gff
    path "${tag}.${outname}.merged.${params.freq}.ffn", emit: mpd_ffn
    path "${tag}.${outname}.merged.${params.freq}.faa", emit: mpd_faa
    path "${tag}.${outname}.merged.${params.freq}.genes", emit: mpd_out
  script:
  """
  prodigal -i ${shortFile.baseName}.fq -d ${tag}.${outname}.short.${params.freq}.ffn -f gff -o ${tag}.${outname}.short.${params.freq}.gff -p meta -q -s ${tag}.${outname}.short.${params.freq}.genes -a ${tag}.${outname}.short.${params.freq}.faa
  prodigal -i ${longFile.baseName}.fq -d ${tag}.${outname}.long.${params.freq}.ffn -f gff -o ${tag}.${outname}.long.${params.freq}.gff -p single -q -s ${tag}.${outname}.long.${params.freq}.genes -a ${tag}.${outname}.long.${params.freq}.faa
  cat ${tag}.${outname}.long.${params.freq}.genes ${tag}.${outname}.short.${params.freq}.genes > ${tag}.${outname}.merged.${params.freq}.genes
  cat ${tag}.${outname}.long.${params.freq}.gff ${tag}.${outname}.short.${params.freq}.gff > ${tag}.${outname}.merged.${params.freq}.gff
  cat ${tag}.${outname}.long.${params.freq}.ffn ${tag}.${outname}.short.${params.freq}.ffn > ${tag}.${outname}.merged.${params.freq}.ffn
  cat ${tag}.${outname}.long.${params.freq}.faa ${tag}.${outname}.short.${params.freq}.faa > ${tag}.${outname}.merged.${params.freq}.faa
  """
}

process callIMPP{
  //publishDir "${params.tmpdir}"
  publishDir "${params.outdir}", mode: 'copy'
  input:
  path(mergedFile)
  path(edge_gff)
  val(outname)
  val(tag)

  output:
  path "${tag}.${outname}.fa", emit: impp_paths

  script:
  """
  impp_extendAnchor -g $mergedFile -a $edge_gff -l $params.maxlen  -p ${tag}.${outname}.fa -k $params.kmin
  """
}

process callBWA{
  executor 'local'
  cpus params.threads
  cache 'deep'
  //publishDir "${params.tmpdir}"
  publishDir "${params.outdir}", mode: 'copy'
  
  input:
  path(mergedFile)
  path(imppFile)
  path(readsfq)
  val(paired)
  val(tag)

  output:
  path "${tag}.edges.${params.bwa_min_score}.sam", emit: bwa_edge
  path "${tag}.paths.${params.bwa_min_score}.sam", emit : bwa_path

  script:
  if(paired == 1)
  """
  bwa index -p edges.index $mergedFile
  bwa mem edges.index -p $readsfq -T $params.bwa_min_score -t $params.threads -k $params.bwa_min_seed -o ${tag}.edges.${params.bwa_min_score}.sam
  
  bwa index -p paths.index $imppFile
  bwa mem paths.index -p $readsfq -T $params.bwa_min_score -t $params.threads -k $params.bwa_min_seed -o ${tag}.paths.${params.bwa_min_score}.sam
  """
  else
  """
  bwa index -p edges.index $mergedFile
  bwa mem edges.index $readsfq -T $params.bwa_min_score -t $params.threads -k $params.bwa_min_seed -o ${tag}.edges.${params.bwa_min_score}.sam
  
  bwa index -p paths.index $imppFile
  bwa mem paths.index $readsfq -T $params.bwa_min_score -t $params.threads -k $params.bwa_min_seed -o ${tag}.paths.${params.bwa_min_score}.sam
  """

}

process sixFrameTranslate{
  cache 'deep'
  //publishDir "${params.tmpdir}"
  publishDir "${params.outdir}", mode: 'copy'
  container 'staphb/spades:latest'
  
  input:
  path(reads_fa) 
  path(reads_gff)
  path(edges_gff)
  path(paths_gff)
  path(edges_bwa)
  path(paths_bwa)
  val(tag)

  output:
  path "${tag}.reads.filter.${params.bwa_min_score}.fasta", emit: filtered_reads
  path "${tag}.reads.translated.${params.bwa_min_score}.fasta", emit : translated_reads
  script:
  """
  sixFrameTranslate.py $reads_gff $edges_gff $paths_gff $edges_bwa $paths_bwa $reads_fa $params.outdir
  mv reads.filter.fasta ${tag}.reads.filter.${params.bwa_min_score}.fasta
  mv reads.translated.fasta ${tag}.reads.translated.${params.bwa_min_score}.fasta
  """

}

/** first round peptide assembly **/
process peptideAssembly_first{
  executor 'local'
  cpus params.threads
   
  cache 'deep'
  //publishDir "${params.tmpdir}"
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(filtered_reads)
  val(tag)

  output:
  path "${tag}.assembled_proteins.${params.bwa_min_score}.faa", emit: plass_first_round

  script:
  """
  plass assemble --threads $params.threads --min-length $params.plass_min_length --num-iterations $params.plass_num_iter $filtered_reads ${tag}.assembled_proteins.${params.bwa_min_score}.faa plass.tmp
  rm -rf plass.tmp
  """
}

process modifyAssemblyOut{
  executor 'local'
  cpus params.threads
  cache 'deep'
  container 'staphb/spades:latest'
  
  //publishDir "${params.tmpdir}"
  publishDir "${params.outdir}", mode: 'copy'
  input:
  path(assembled_proteins)
  val(tag)
  output:
  path "${tag}.assembled_proteins.${params.bwa_min_score}.re.60.faa", emit: plass_renamed
  """
  renamePlassContig.py $assembled_proteins
  """

}
process alignSixFrametoAssembly{
  executor 'local'
  cpus params.threads
  
  cache 'deep'
  //publishDir "${params.tmpdir}"
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(assembled_proteins)
  path(translated_reads)
  val(tag)
  output:
  path "${tag}.plass.translated.${params.bwa_min_score}.dmd", emit: plass_alignment
  script:
  """
  diamond makedb --in $assembled_proteins -d ${tag}.plass_contig.index
  diamond blastx -p $params.threads -d ${tag}.plass_contig.index -q $translated_reads -o ${tag}.plass.translated.${params.bwa_min_score}.dmd
  """
}

process refinePlassInput{
  executor 'local'
  cpus params.threads  
  cache 'deep'
  container 'staphb/spades:latest'
  //publishDir "${params.tmpdir}"
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(assembled_proteins)
  path(plass_dmd)
  path(translated_reads)
  path(filtered_reads)
  val(tag)
  output:
  path "${tag}.reads.impp.fasta", emit: impp_filtered_reads
  path "${tag}.reads.translated.${params.bwa_min_score}.mapped.fasta", emit: impp_mapped_reads

  script:
  """
  plassMapping.py $assembled_proteins $plass_dmd $translated_reads  
  getPredReads.py $filtered_reads $translated_reads  
  cat ${tag}.reads.filter.${params.bwa_min_score}.pred.fasta ${tag}.reads.translated.${params.bwa_min_score}.mapped.fasta > ${tag}.reads.impp.fasta
  """
}

/** second round of plass assembly **/
process peptideAssembly_second{
  executor 'local'
  cpus params.threads 
  cache 'deep'
  
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(impp_reads)
  val(tag)

  output:
  path "${tag}.assembled_proteins.${params.bwa_min_score}.impp.faa", emit: final_plass_assembly

  script:
  """
  plass assemble --threads $params.threads --min-length $params.plass_min_length --num-iterations $params.plass_num_iter $impp_reads ${tag}.assembled_proteins.${params.bwa_min_score}.impp.faa plass.tmp
  rm -rf plass.tmp
  """
}





   


