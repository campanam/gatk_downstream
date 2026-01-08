#!/usr/bin/env nextflow

/* GATK Joint-Genotyping Pipeline version 0.2.0
Michael G. Campana, 2023-2026
Smithsonian\'s National Zoo and Conservation Biology Institute

The software is made available under the Smithsonian Institution terms of use (https://www.si.edu/termsofuse). */

gatk = 'gatk --java-options "' + params.java_options + '" ' // Simplify gatk command line

process refDictFai {
	
	// Prepare nuclear reference sequence dict and fai indices
	
	input:
	path refseq
	
	output:
	path "${refseq.baseName}*.{fai,dict}"

	"""
	samtools faidx ${refseq}
	samtools dict ${refseq} > ${refseq.baseName}.dict
	"""

}

process extractChrNames {

	// Extract chromosome names from reference sequence to use all chromosomes
	
	input:
	path refseq
	
	output:
	path "${refseq.baseName}.chr.txt"
	
	"""
	grep '>' $refseq | sed 's/>//g' | cut -f1 -d ' ' > ${refseq.baseName}.chr.txt
	"""

}

process createGenomicsDB {

	// Create chromosome-specific GenomicsDB
	
	publishDir "$params.outdir/01_ChrGenomicsDBs", mode: 'copy'
	
	input:
	path gvcfs
	val chr
	val stem
	
	output:
	path "${stem}_${chr}.tgz"
	
	"""
	VARPATH=""
	for file in ${gvcfs}/*.vcf.gz; do VARPATH+=" -V \$file"; done
	$gatk GenomicsDBImport\$VARPATH --genomicsdb-workspace-path ${stem}_${chr} --intervals $chr
	tar czf ${stem}_${chr}.tgz ${stem}_${chr}/*
	"""

}

process genMapIndex {

	// Generate GenMap index. From RatesTools 1.2.4
	
	label 'genmap'
		
	input:
	path refseq
	val gm_tmpdir
	
	output:
	tuple path("$refseq"), path("${refseq.simpleName}_index"), path("${refseq.simpleName}_index/*")
	
	"""
	export TMPDIR=${gm_tmpdir}
	if [ ! -d ${gm_tmpdir} ]; then mkdir ${gm_tmpdir}; fi
	genmap index -F ${refseq} -I ${refseq.simpleName}_index
	"""

}

process genMapMap {

	// Calculate mappability using GenMap and filter using filterGM. From RatesTools 1.2.4
	
	label 'genmap'
	label 'ruby'
		
	input:
	tuple path(refseq), path(genmap_index), path("*")
	
	output:
	path "${refseq.simpleName}_genmap.1.0.bed"
	
	"""
	genmap map ${params.gm_opts} -T ${task.cpus} -I ${refseq.simpleName}_index/ -O ${refseq.simpleName}_genmap -b
	filterGM.rb ${refseq.simpleName}_genmap.bed 1.0 exclude > ${refseq.simpleName}_genmap.1.0.bed
	"""
}

workflow.onComplete {
	if (workflow.success) {
		println "Joint-genotyping pipeline completed successfully at $workflow.complete!"
		if (params.email != "NULL") {
			sendMail(to: params.email, subject: 'Joint-genotyping pipeline successful completion', body: "Joint-genotyping pipeline completed successfully at $workflow.complete!")
		}
	} else {
		println "Joint-genotyping pipeline terminated with errors at $workflow.complete.\nError message: $workflow.errorMessage"
		if (params.email != "NULL") {
			sendMail(to: params.email, subject: 'Joint-genotyping pipeline terminated with errors', body: "Joint-genotyping pipeline terminated with errors at $workflow.complete.\nError message: $workflow.errorMessage")
		}
	}
}

workflow {
	main:
		refFaiDict(params.refseq)

}

	