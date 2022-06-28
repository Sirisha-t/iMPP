#! /usr/bin/env nextflow

nextflow.enable.dsl=2

/** 
Input parameters
**/
params.interleaved = null
params.forward = null
params.reverse = null
params.single = null
params.outdir = "${workflow.projectDir}/output"
params.genecaller = "fgs"
params.threads = 8
params.maxlen = 150
params.help = false

/** Print help message **/
def helpMessage(){
log.info """
 USAGE: nextflow run main.nf --single [other options] <fastq file/s> --outdir <output_directory> -profile docker
 Input data options:
   -profile               <string>      : [required] docker 
   --interleaved          <filename>    : fastq file with interlaced forward and reverse paired-end reads
   --forward              <filename>    : fastq file with forward paired-end reads
   --reverse              <filename>    : fastq file with reverse paired-end reads
   --single               <filename>    : fastq file with unpaired reads
   --outdir               <dirname>     : [required] output directory name including the full path [default: output]
   --genecaller           <string>      : [optional] fgs or prodigal [default: fgs]
   --threads              <int>         : [optional] number of threads [default: 16]
   --maxlen               <int>         : [optional] maximum extension length for anchors [default: 150]
   -h/--help                            :  help message

  NOTE: Input FASTQ read file must be specified in either --interleaved, --forward & --reverse or --single format.

"""
}

// Show help message
if (params.help){
  helpMessage()
  exit 0
}

workflow
{
  /** Checking if input parameters have been entered correctly **/
  if (params.single == "null" && params.forward=="null" && params.reverse=="null" && params.interleaved=="null"){
    println("\nERROR: Input fastq files were not specified properly. Please see usage below. \n\n")
    helpMessage()
    exit 1
  }

  outpath = file(params.outdir)
  if(!outpath.exists())
  {
    result = outpath.mkdir()
    println result ? "Created output dir" : "Cannot create directory: $outpath"
  }

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
            if( (params.forward.getExtension() ==  "fq" || params.forward.getExtension() ==  "fastq") && (params.reverse.getExtension() ==  "fq" || params.reverse.getExtension() ==  "fastq")){
              fwd_file = Channel.fromPath(params.forward)
              rev_file = Channel.fromPath(params.reverse)
              mergeFastqFiles(fwd_file, rev_file)
              inputfile = Channel.fromPath(mergedFastqFiles.out.merged_file)
              outputfile = Channel.fromPath(outpath)
              convertFastq(inputfile, outputfile) 
              pairend = 0
            }
            else{
              println("\nERROR: Input files must end in .fq or .fastq extension. Please unzip any zipped input files and try again.\n\n")
              exit 1
            }
          }

        }
      else{
          if( (params.interleaved.getExtension() ==  "fq"|| params.interleaved.getExtension() == "fastq")){
            fastq_file = file(params.interleaved)
            inputfile = Channel.fromPath(params.interleaved)
            outputfile = Channel.fromPath(outpath)
            convertFastq(inputfile, outputfile) 
            convertFastq()
            paired = 0
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
      }
      else { 
          println("\nERROR: Input files must end in .fq or .fastq extension. Please unzip any zipped input files and try again.\n\n")
          exit 1
      }
    }

  if (params.maxlen && params.maxlen >= 350 ){
      println("\nERROR: Maximum extend length threshold too high. Please select a lower value for optimal results. [default=150] \n\n")
      helpMessage()
      exit 1
    }

  if(params.genecaller == "fgs")
    {
    callFGS(convertFastq.out.reads_fasta)
    genecaller = 1
    }

  else if(params.genecaller == "prodigal")
    {
    callMPD(convertFastq.out.reads_fasta)
    genecaller = 0
    }
    
  else{
    println("ERROR: Invalid gene caller option specified. Please use either fgs or prodigal. [default=fgs]")
    helpMessage()
    exit 1
    }

    callSGAAssembly(convertFastq.out.reads_fastq)
    callSPAdesAssembly(convertFastq.out.reads_fastq)
    getMergedGraph(callSGAAssembly.out.asqg_file, convertFastq.out.reads_fastq, callFGS.out.reads_gff, callSPAdesAssembly.out.contig)

    if(genecaller){
    edge_out = "edges"
    splitEdges(getMergedGraph.out.merged_graph)
    callFGSEdges(splitEdges.out.short_edge, splitEdges.out.long_edge, edge_out)
    
    path_out = "paths"
    callIMPP(getMergedGraph.out.merged_graph, callFGSEdges.out.fgs_edge_gff, path_out)
    splitPaths(callIMPP.out.impp_paths)
    callFGSPaths(splitPaths.out.short_path, splitPaths.out.long_path, path_out)
    }
    else{
    }

    callBWA(getMergedGraph.out.merged_graph, callIMPP.out.impp_paths, convertFastq.out.reads_fastq)
    
    /** Last step -- perform 2 rounds of peptide assembly **/
    sixFrameTranslate(convertFastq.out.reads_fasta, callFGS.out.reads_gff, callFGSEdges.out.fgs_edge_gff, callFGSPaths.out.fgs_path_gff, callBWA.out.bwa_edge, callBWA.out.bwa_path)
    peptideAssembly_R1(sixFrameTranslate.out.filtered_reads)
    modifyAssemblyOut(peptideAssembly_R1.out.plass_first_round)
    alignSixFrametoAssembly(modifyAssemblyOut.out.plass_renamed, sixFrameTranslate.out.translated_reads)
    refinePlassInput(modifyAssemblyOut.out.plass_renamed, alignSixFrametoAssembly.out.plass_alignment, sixFrameTranslate.out.translated_reads, sixFrameTranslate.out.filtered_reads )
    peptideAssembly_R2(refinePlassInput.out.impp_filtered_reads)
    
    cleanOutput(callFGS.out.reads_gff, callFGS.out.reads_ffn, callFGS.out.reads_faa, callFGSPaths.out.fgs_path_gff, callFGSPaths.out.fgs_path_ffn, callFGSPaths.out.fgs_path_faa, peptideAssembly_R2.out.final_plass_assembly)    

}


process convertFastq{
  cache 'deep'
  publishDir "${params.outdir}"
  input:
  path queryFile
  path out

  output:
  path 'reads.fasta', emit: reads_fasta
  path 'reads.fastq', emit: reads_fastq
  """
  seqtk seq -A $queryFile > query.fasta
  seqtk rename query.fasta > reads.fasta
  seqtk rename $queryFile > reads.fastq
  """
}


process mergeFastqFiles{
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

process callFGS{
  cache 'deep'
  executor 'local'
  cpus params.threads
  memory '40 GB'
  input:
    path(queryFile)
  output:
    publishDir "${params.outdir}"
    path 'reads.gff', emit: reads_gff
    path 'reads.ffn', emit: reads_ffn
    path 'reads.faa', emit: reads_faa

  script:
  """
  run_FragGeneScan.pl -thread=$params.threads -genome=$queryFile -out=reads -complete=$params.fgs_mode -train=$params.fgs_train
  """
}


process callMPD{
  executor 'local'
  cpus params.threads
  memory '40 GB'  
  input:
    path(queryFile)
  output:
    publishDir "${params.outdir}"
    path 'reads.*', emit: mpd_reads
  script:
  """
  prodigal -i $queryFile -d reads.ffn -f gff -o reads.gff -p meta -q -s reads.genes -a reads.faa
  """
}


process callSGAAssembly{
  cache 'deep'
  executor 'local'
  cpus params.threads
  memory '40 GB'  
  input:
  path(queryFile)
  output:
  path 'assembly.*', emit: sga_assembly
  path 'assembly.ec.rmdup.asqg', emit: asqg_file
  publishDir "${params.outdir}"
  script:
  """
  sga preprocess -o assembly.pp.fq $queryFile
  sga index -t $params.threads assembly.pp.fq
  sga correct -t $params.threads -k $params.sga_correct_kmer --learn -m $params.sga_correct_min_olap assembly.pp.fq -o assembly.ec.fq
  sga index -t $params.threads -a ropebwt assembly.ec.fq
  sga rmdup -t $params.threads assembly.ec.fq
  sga overlap -t $params.threads -m $params.sga_overlap_min_len assembly.ec.rmdup.fa
  gunzip -f assembly.ec.rmdup.asqg.gz
  """
}

process getMergedGraph{
  cache 'deep'
  container 'biocontainers/bwa:v0.7.17_cv1'
  publishDir "${params.outdir}"
  input:
  path(asqg)
  path(queryFile)
  path(reads_gff)
  path(spades_contig)
  output:
  path 'og.merged.fq', emit: merged_graph
  """
  impp_modASQG -g $asqg
  impp_getSG -g assembly.ec.rmdup.mod.asqg -r $queryFile -f $reads_gff
  mv assembly.ec.rmdup.mod.StringGraph.fq  og.fq

  bwa index -p contigs.index $spades_contig
  bwa mem contigs.index og.fq -T 45 -t $params.threads -k $params.bwa_min_seed -o og.sam

  impp_mergeSG og.fq og.sam $spades_contig
  """
} 


process callSPAdesAssembly{
  cache 'deep'
  executor 'local'
  cpus params.threads
  memory '40 GB' 
  publishDir "${params.outdir}"
  input:
  path(queryFile)
  output:
  path 'spades/assembly_graph.dbGraph.*', emit: spades_assembly
  path 'spades/contigs.fasta', emit: contig

  script:
  """
  spades.py -t $params.threads -m $params.spades_memory_limit -o $params.spadesdir -s $queryFile --only-assembler -k $params.spades_kmer 
  renameDB.py $params.spadesdir/assembly_graph.fastg
  """
}

process splitEdges{
  cache 'deep'
  container 'staphb/spades:latest'
  publishDir "${params.outdir}"
  input:
  path(mergedFile)
  output:
  path "${mergedFile.baseName}.short.*", emit: short_edge
  path "${mergedFile.baseName}.long.*", emit: long_edge

  """
  splitFiles.py $mergedFile $params.freq 
  """

}
process splitPaths{
  cache 'deep'
  container 'staphb/spades:latest'
  publishDir "${params.outdir}"
  input:
  path(mergedFile)
  output:
  path "${mergedFile.baseName}.short.*", emit: short_path
  path "${mergedFile.baseName}.long.*", emit: long_path

  """
  splitFiles.py $mergedFile $params.freq 
  """

}

process callFGSEdges{
  executor 'local'
  cpus params.threads
  memory '40 GB' 
  container 'tsirisha/fraggenescan:latest'

  publishDir "${params.outdir}"
  input:
  path(shortFile)
  path(longFile)
  val(outname)
  output:
    path "${outname}.merged.${params.freq}.gff", emit: fgs_edge_gff
    path "${outname}.merged.${params.freq}.ffn", emit: fgs_edge_ffn
    path "${outname}.merged.${params.freq}.faa", emit: fgs_edge_faa
  script:
  """
  
  run_FragGeneScan.pl -genome=${shortFile.baseName}.fq -out=${outname}.short.${params.freq} -complete=0 -train=$params.fgs_train_other -thread=$params.threads
  run_FragGeneScan.pl -genome=${longFile.baseName}.fq -out=${outname}.long.${params.freq} -complete=1 -train=$params.fgs_train_complete -thread=$params.threads
  cat ${outname}.long.${params.freq}.out ${outname}.short.${params.freq}.out > ${outname}.merged.${params.freq}.out
  cat ${outname}.long.${params.freq}.gff ${outname}.short.${params.freq}.gff > ${outname}.merged.${params.freq}.gff
  cat ${outname}.long.${params.freq}.ffn ${outname}.short.${params.freq}.ffn > ${outname}.merged.${params.freq}.ffn
  cat ${outname}.long.${params.freq}.faa ${outname}.short.${params.freq}.faa > ${outname}.merged.${params.freq}.faa
  """
}

process callFGSPaths{
  executor 'local'
  cpus params.threads
  memory '40 GB' 
  container 'tsirisha/fraggenescan:latest'

  publishDir "${params.outdir}"
  input:
  path(shortFile)
  path(longFile)
  val(outname)
  output:
  path "${outname}.merged.${params.freq}.gff", emit: fgs_path_gff
  path "${outname}.merged.${params.freq}.ffn", emit: fgs_path_ffn
  path "${outname}.merged.${params.freq}.faa", emit: fgs_path_faa

  script:
  """
  
  run_FragGeneScan.pl -genome=${shortFile.baseName}.fq -out=${outname}.short.${params.freq} -complete=0 -train=$params.fgs_train_other -thread=$params.threads
  run_FragGeneScan.pl -genome=${longFile.baseName}.fq -out=${outname}.long.${params.freq} -complete=1 -train=$params.fgs_train_complete -thread=$params.threads
  cat ${outname}.long.${params.freq}.out ${outname}.short.${params.freq}.out > ${outname}.merged.${params.freq}.out
  cat ${outname}.long.${params.freq}.gff ${outname}.short.${params.freq}.gff > ${outname}.merged.${params.freq}.gff
  cat ${outname}.long.${params.freq}.ffn ${outname}.short.${params.freq}.ffn > ${outname}.merged.${params.freq}.ffn
  cat ${outname}.long.${params.freq}.faa ${outname}.short.${params.freq}.faa > ${outname}.merged.${params.freq}.faa
  """
}

process callIMPP{
  publishDir "${params.outdir}"
  input:
  path(mergedFile)
  path(edge_gff)
  val(outname)

  output:
  path "${outname}.fa", emit: impp_paths

  script:
  """
  impp_extendAnchor -g $mergedFile -a $edge_gff -l $params.maxlen  -p ${outname}.fa -k $params.kmin
  """
}

process callBWA{
  executor 'local'
  cpus params.threads
  memory '40 GB' 
  cache 'deep'
  publishDir "${params.outdir}"
  input:
  path(mergedFile)
  path(imppFile)
  path(readsfq)

  output:
  path "edges.${params.bwa_min_score}.sam", emit: bwa_edge
  path "paths.${params.bwa_min_score}.sam", emit : bwa_path

  script:
  """
  bwa index -p edges.index $mergedFile
  bwa mem edges.index $readsfq -T $params.bwa_min_score -t $params.threads -k $params.bwa_min_seed -o edges.${params.bwa_min_score}.sam
  
  bwa index -p paths.index $imppFile
  bwa mem paths.index $readsfq -T $params.bwa_min_score -t $params.threads -k $params.bwa_min_seed -o paths.${params.bwa_min_score}.sam
  """

}

process sixFrameTranslate{
  cache 'deep'
  publishDir "${params.outdir}"
  input:
  path(reads_fa) 
  path(reads_gff)
  path(edges_gff)
  path(paths_gff)
  path(edges_bwa)
  path(paths_bwa)

  output:
  path "reads.filter.${params.bwa_min_score}.fasta", emit: filtered_reads
  path "reads.translated.${params.bwa_min_score}.fasta", emit : translated_reads
  script:
  """
  python3 $projectDir/bin/sixFrameTranslate.py $reads_gff $edges_gff $paths_gff $edges_bwa $paths_bwa $reads_fa $params.outdir
  mv reads.filter.fasta reads.filter.${params.bwa_min_score}.fasta
  mv reads.translated.fasta reads.translated.${params.bwa_min_score}.fasta
  """

}

/** first round peptide assembly **/
process peptideAssembly_R1{
  executor 'local'
  cpus params.threads
  memory '40 GB' 
  cache 'deep'
  publishDir "${params.outdir}"

  input:
  path(filtered_reads)

  output:
  path "assembled_proteins.${params.bwa_min_score}.faa", emit: plass_first_round

  script:
  """
  plass assemble --threads $params.threads --min-length $params.plass_min_length --num-iterations $params.plass_num_iter $filtered_reads assembled_proteins.${params.bwa_min_score}.faa plass.tmp
  rm -rf plass.tmp
  """
}

process modifyAssemblyOut{
  executor 'local'
  cpus params.threads
  memory '40 GB' 
  cache 'deep'
  publishDir "${params.outdir}"
  input:
  path(assembled_proteins)
  output:
  path "assembled_proteins.${params.bwa_min_score}.re.60.faa", emit: plass_renamed
  """
  python $projectDir/bin/renamePlassContig.py $assembled_proteins
  """

}
process alignSixFrametoAssembly{
  executor 'local'
  cpus params.threads
  memory '40 GB' 
  cache 'deep'
  publishDir "${params.outdir}"

  input:
  path(assembled_proteins)
  path(translated_reads)
  output:
  path "plass.translated.${params.bwa_min_score}.dmd", emit: plass_alignment
  script:
  """
  diamond makedb --in $assembled_proteins -d plass_contig.index
  diamond blastx -p $params.threads -d plass_contig.index -q $translated_reads -o plass.translated.${params.bwa_min_score}.dmd
  """
}

process refinePlassInput{
  executor 'local'
  cpus params.threads
  memory '40 GB' 
  cache 'deep'
  publishDir "${params.outdir}"

  input:
  path(assembled_proteins)
  path(plass_dmd)
  path(translated_reads)
  path(filtered_reads)
  output:
  path 'reads.impp.fasta', emit: impp_filtered_reads
  script:
  """
  python $projectDir/bin/plassMapping.py $assembled_proteins $plass_dmd $translated_reads  
  python $projectDir/bin/getPredReads.py $filtered_reads $translated_reads  
  cat reads.filter.${params.bwa_min_score}.pred.fasta reads.translated.${params.bwa_min_score}.mapped.fasta > reads.impp.fasta
  """
}

/** second round of plass assembly **/
process peptideAssembly_R2{
  executor 'local'
  cpus params.threads
  memory '40 GB' 
  cache 'deep'
  publishDir "${params.outdir}"

  input:
  path(impp_reads)

  output:
  path "assembled_proteins.${params.bwa_min_score}.impp.faa", emit: final_plass_assembly

  script:
  """
  plass assemble --threads $params.threads --min-length $params.plass_min_length --num-iterations $params.plass_num_iter $impp_reads assembled_proteins.${params.bwa_min_score}.impp.faa plass.tmp
  rm -rf plass.tmp
  """
}

process cleanOutput{
  publishDir "${params.outdir}", mode: 'move'
  input:
  path(reads_gff)
  path(reads_ffn)
  path(reads_faa)
  path(paths_gff)
  path(paths_ffn)
  path(paths_faa)
  path(assembled_proteins)

  output:
  path 'orfs.*', emit: final_output
  """
  cat $reads_gff $paths_gff > orfs.gff
  cat $reads_ffn $paths_ffn > orfs.ffn
  cat $reads_faa $paths_faa > orfs.faa
  cat $assembled_proteins > assembled_proteins.faa
  """
}

workflow.onComplete {
    println "iMPP Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "ERROR: iMPP Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

