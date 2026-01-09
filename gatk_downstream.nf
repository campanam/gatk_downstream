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

process jointGenotype {

	// Joint-genotype gVCFs
	
	publishDir "$params.outdir/02_RawChrVCFs", mode: 'copy'
	
	input:
	path db
	path refseq
	path "*"
	
	output:
	path "${db.baseName}.vcf.gz"
	
	"""
	tar xfz $db
	$gatk GenotypeGVCFs -R $refseq -V gendb://${db.baseName} -O ${db.baseName}.vcf
	bgzip ${db.baseName}.vcf
	rm -r ${db.baseName}
	"""

}

process vcftoolsSiteFilter {

	// Perform VCFtools site filters on VCFs
	// Modified from RatesTools 1.2.4: Armstrong & Campana 2023
	
	publishDir "$params.outdir/03_VCFtoolsFilterChrVCFs", mode: 'copy', pattern: '*vcftools.vcf.gz'
	
	input:
	path raw_vcf
	
	output:
	path("${raw_vcf.baseName[0..-5]}.vcftools.tmp")
	path(raw_vcf)
	path("${raw_vcf.baseName[0..-5]}.vcftools.vcf.gz")
	
	script:
	if (params.vcftools_site_filters == "NULL")
		"""
		cp -P $raw_vcf ${raw_vcf.baseName[0..-5]}.vcftools.vcf.gz
		vcftools --gzvcf $raw_vcf
		cp .command.log ${raw_vcf.baseName[0..-5]}.vcftools.tmp
		"""
	else
		"""
		vcftools --gzvcf ${raw_vcf} --recode -c ${params.vcftools_site_filters} | bgzip > ${raw_vcf.baseName[0..-5]}.vcftools.vcf.gz
		cp .command.log ${raw_vcf.baseName[0..-5]}.vcftools.tmp
		"""

}

process sanityCheckLogs {

	// Sanity check filtering logs and remove too short contigs as needed

	label 'gzip'
	
	input:
	path logfile
	path allvcflog
	path filtvcflog
	val min_contig_length
	val min_filt_contig_length
	
	output:
	path "${logfile.baseName}.log",  emit: log
	path "${filtvcflog.baseName.split(".vcf")[0]}.OK.vcf.gz", optional: true, emit: ok_vcf
	
	"""
	logstats.sh $logfile $allvcflog $filtvcflog $min_contig_length $min_filt_contig_length  > ${logfile.baseName}.log
	"""
	
}
	
	
process gatkSiteFilter {

	// Perform GATK site filters on VCFs
	// Modified from RatesTools 1.2.4: Armstrong & Campana 2023
	
	publishDir "$params.outdir/04_GATKFilterChrVCFs", mode: 'copy', pattern: '*gatk.vcf.gz'
	
	input:
	path vcftools_vcf
	path refseq from params.refseq
	path "*"
	
	output:
	path("${vcftools_vcf.baseName[0..-17]}.gatk.tmp")
	path(vcftools_vcf)
	path("${vcftools_vcf.baseName[0..-17]}.gatk.vcf.gz")
	
	script:
	if (params.gatk_site_filters == "NULL")
		"""
		ln -s $vcftools_vcf ${vcftools_vcf.baseName[0..-17]}.gatk.vcf.gz
		vcftools --gzvcf $vcftools_vcf
		cp .command.log ${vcftools_vcf.baseName[0..-17]}.gatk.tmp
		"""
	else
		"""
		tabix $vcftools_vcf
		$gatk VariantFiltration -R $refseq -V $vcftools_vcf -O tmp.vcf.gz ${params.gatk_site_filters}
		$gatk SelectVariants -R $refseq -V tmp.vcf.gz -O ${vcftools_vcf.baseName[0..-17]}.gatk.vcf.gz --exclude-filtered
		rm tmp.vcf.gz
		vcftools --gzvcf ${vcftools_vcf.baseName[0..-17]}.gatk.vcf.gz
		tail .command.log > ${vcftools_vcf.baseName[0..-17]}.gatk.tmp
		"""

}

process concatenateVCFs {

	// Concatenate VCFs using BCFtools concat
	// Attempt to sort by order in original file
	
	publishDir "$params.outdir/05_ConcatVCF", mode: 'copy'
	
	input:
	path chrfile
	path "*"
	
	output:
	path "${params.stem}.all.vcf.gz"
	
	"""
	#!/usr/bin/env bash
	fileline=''
	readarray -t chrs < $chrfile
	len=\${#chrs[@]}
	for ((i=0; i<\$len; i++)); do
		if [ -f ${params.stem}_\${chrs[\$i]}.gatk.OK.vcf.gz ]; then
			fileline+="  ${params.stem}_\${chrs[\$i]}.gatk.OK.vcf.gz"
		fi
	done
	bcftools concat -O v -o ${params.stem}.all.vcf.gz\$fileline
	"""

}

process filterMappability {

	// Remove sites with mappability < 1.0 using BEDtools
	
	publishDir "$params.outdir/06_MapFiltVCF", mode: 'copy'
	
	input:
	path vcf
	path bed
	
	output:
	path "${params.stem}.map.vcf.gz"
	
	"""
	bedtools subtract -a $vcf -b $bed -header | gzip > ${params.stem}.map.vcf.gz
	"""

}

process snpRelate {

	// Estimate kinship using SNPRelate
	
	publishDir "$params.outdir/07_SNPRelate", mode: 'copy'
	
	input:
	path vcf from map_vcf_ch
	
	output:
	path "${params.stem}.gds"
	path "${params.stem}*.csv"
	path "${params.stem}.snprelate.log"
	path "${params.stem}.Rdata"
	
	"""
	#!/usr/bin/env Rscript
	library("SNPRelate")
	source(system("which kinshipUtils.R", intern = TRUE))
	snpgdsVCF2GDS(Sys.readlink(\'$vcf\'), \'${params.stem}.gds\', method = "biallelic.only")
	snps <- snpgdsOpen(\'${params.stem}.gds\')
	pruned <- snpgdsLDpruning(snps, ${params.snprelate_opts}, ld.threshold = ${params.snprelate_ld})
	whole_kinship <- snpgdsIBDMLE(snps, snp.id = unlist(pruned), ${params.snprelate_opts}, num.thread = ${task.cpus})
	bootstrapped <- bootstrap.kinship(snps, ibdmethod = "MLE", ${params.snprelate_opts}, num.thread = ${task.cpus}, ld.threshold = ${params.snprelate_ld})
	write.kinship.matrix(bootstrapped, meanfile = \"${params.stem}_bootstrap_meanvalues.csv\", cifile = \"${params.stem}_random_kinship_CI.csv\")
	system(\"cp .command.log ${params.stem}.snprelate.log\")
	save.image(file = \"${params.stem}.Rdata\")
	"""

}

process ngsRelate {

	// Estimate kinship using ngsRelate
	// Requires local compilation and installation of ngsRelate v2 commit ec95c8f built with htslib v 1.18 commit 99415e2
	
	publishDir "$params.outdir/08_ngsRelate", mode: 'copy'
	
	input:
	path vcf
	
	output:
	path "${params.stem}.ngsRelate.res"
	
	"""
	zgrep -m1 '#CHROM' $vcf | cut -f10- | sed "s/\\t/\\n/g" > sampleids.txt
	ngsRelate -h $vcf -O ${params.stem}.ngsRelate.res ${params.ngsrelate_opts} -p ${task.cpus} -z sampleids.txt -n `wc -l sampleids.txt`
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

workflow logVcftoolsSanity {
	// Sanity check logs from VCFtools site filtering. From RatesTools 1.2.4
	take:
		tmpfile
		rawvcf
		filtvcf
	main:
		sanityCheckLogs(tmpfile, rawvcf, filtvcf, params.min_contig_length, params.min_filt_contig_length)
	emit:
		sanelog = sanityCheckLogs.out.log
		ok_vcf = sanityCheckLogs.out.ok_vcf
}

workflow logGatkSanity {
	// Sanity check logs for GATK site filtering and remove too short contigs. From RatesTools 1.2.4
	// Dummy value of 1 for min_contig_length since already evalutated and no longer accurate
	take:
		tmpfile
		rawvcf
		filtvcf
	main:
		sanityCheckLogs(tmpfile, rawvcf, filtvcf, 1, params.min_filt_contig_length)
	emit:
		sanelog = sanityCheckLogs.out.log
		ok_vcf = sanityCheckLogs.out.ok_vcf
}


workflow {
	main:
		refDictFai(params.refseq)
		genMapIndex(params.refseq, params.gm_tmpdir) | genMapMap
		if (params.chrlist == "NULL") {
			extractChrNames(params.refseq)
			chr_ch = extractChrNames.out.flatten().splitText(by: 1).map { it.replaceAll(/\n/, "") }
		} else {
			chr_ch = channel.fromPath(params.chrlist).flatten().splitText(by: 1).map { it.replaceAll(/\n/, "") }
		}
		createGenomicsDB(params.gvcfs, chr_ch, params.stem)
		jointGenotype(createGenomicsDB.out, params.refseq, refDictFai.out) | vcftoolsSiteFilter | logVcftoolsSanity
		gatkSiteFilter(logVcftoolsSanity.out.ok_vcf, params.refseq, refDictFai.out) | logGatkSanity
		if (params.chrlist == "NULL") {
			concatenateVCFs(extractChrNames.out, logGatkSanity.out.ok_vcf.collect())
		} else {
			concatenateVCFs(params.chrlist, logGatkSanity.out.ok_vcf.collect())
		}
		filterMappability(concatenateVCFs.out, genMapMap.out) | snpRelate
		ngsRelate(filterMappability.out)
}