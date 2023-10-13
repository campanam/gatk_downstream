#!/usr/bin/env nextflow

gatk = 'gatk --java-options "' + params.java_options + '" ' // Simplify gatk command line

nextflow.enable.dsl=1

process extractChrNames {

	// Extract chromosome names from reference sequence
	
	input:
	path refseq from params.refseq
	
	output:
	path "${refseq.baseName}.chr.txt" into chrfile_ch, chrfile_ch2
	
	"""
	grep '>' $refseq | sed 's/>//g' | cut -f1 -d ' ' > ${refseq.baseName}.chr.txt
	"""

}

chr_ch = chrfile_ch.flatten().splitText(by: 1).map { it.replaceAll(/\n/, "") }

process createGenomicsDB {

	// Create chromosome-specific GenomicsDB
	
	publishDir "$params.outdir/01_ChrGenomicsDBs", mode: 'copy'
	
	input:
	path gvcfs from params.gvcfs
	val chr from chr_ch
	val stem from params.stem
	
	output:
	path "${stem}_${chr}.tgz" into genomicsdb_ch
	
	"""
	VARPATH=""
	for file in ${gvcfs}/*.vcf.gz; do VARPATH+=" -V \$file"; done
	$gatk GenomicsDBImport\$VARPATH --genomicsdb-workspace-path ${stem}_${chr} --intervals $chr
	tar czf ${stem}_${chr}.tgz ${stem}_${chr}/*
	"""

}

process genMapIndex {

	// Generate GenMap index. From RatesTools 0.5.16

	
	input:
	path refseq from params.refseq
	val gm_tmpdir from params.gm_tmpdir
	
	output:
	path "${refseq.simpleName}_index" into genmap_index_ch
	path "${refseq.simpleName}_index/*" into genmap_index_files_ch
	
	"""
	export TMPDIR=${gm_tmpdir}
	if [ ! -d ${gm_tmpdir} ]; then mkdir ${gm_tmpdir}; fi
	genmap index -F ${refseq} -I ${refseq.simpleName}_index
	"""

}

process genMapMap {

	// Calculate mappability using GenMap and filter using filterGM.  From RatesTools 0.5.16
	
	publishDir "$params.outdir/06_MapFiltVCF", mode: 'copy'
	
	input:
	path refseq from params.refseq
	path genmap_index from genmap_index_ch
	path '*' from genmap_index_files_ch
	
	output:
	path "${refseq.simpleName}_genmap.1.0.bed" into genmap_ch
	
	"""
	genmap map -K 30 -E 2 -T ${task.cpus} -I ${refseq.simpleName}_index/ -O ${refseq.simpleName}_genmap -b
	filterGM.rb ${refseq.simpleName}_genmap.bed 1.0 exclude > tmp.bed
	bedtools merge -i tmp.bed > ${refseq.simpleName}_genmap.1.0.bed
	"""
}

process jointGenotype {

	// Joint-genotype gVCFs
	
	publishDir "$params.outdir/02_RawChrVCFs", mode: 'copy'
	
	input:
	path db from genomicsdb_ch
	path refseq from params.refseq
	path fai from params.refseq_fai
	path dict from params.refseq_dict
	
	output:
	path "${db.baseName}.vcf.gz" into raw_vcf_ch
	
	"""
	tar xfz $db
	$gatk GenotypeGVCFs -R $refseq -V gendb://${db.baseName} -O ${db.baseName}.vcf
	bgzip ${db.baseName}.vcf
	rm -r ${db.baseName}
	"""

}

process vcftoolsSiteFilter {

	// Perform VCFtools site filters on VCFs
	// Modified from RatesTools 0.5.15: Armstrong & Campana 2023
	
	publishDir "$params.outdir/03_VCFtoolsFilterChrVCFs", mode: 'copy', pattern: '*vcftools.vcf.gz'
	
	input:
	path raw_vcf from raw_vcf_ch
	val site_filters from params.vcftools_site_filters
	
	output:
	tuple path("${raw_vcf.baseName[0..-5]}.vcftools.tmp"), path(raw_vcf), path("${raw_vcf.baseName[0..-5]}.vcftools.vcf.gz")  into vcftools_vcf_ch
	
	script:
	if (site_filters == "NULL")
		"""
		cp -P $raw_vcf ${raw_vcf.baseName[0..-5]}.vcftools.vcf.gz
		vcftools --gzvcf $raw_vcf
		cp .command.log ${raw_vcf.baseName[0..-5]}.vcftools.tmp
		"""
	else
		"""
		vcftools --gzvcf ${raw_vcf} --recode -c ${site_filters} | bgzip > ${raw_vcf.baseName[0..-5]}.vcftools.vcf.gz
		cp .command.log ${raw_vcf.baseName[0..-5]}.vcftools.tmp
		"""

}

process sanityCheckLogsVcftools {

	// Sanity check logs for VCFtools site filtering and remove too short contigs
	// Modified from RatesTools 0.5.15: Armstrong & Campana 2023

	input:
	tuple path(logfile), path(allvcflog), path(filtvcflog) from vcftools_vcf_ch
	val min_contig_length from params.min_contig_length
	val min_filt_contig_length from params.min_filt_contig_length
	
	output:
	path "${logfile.baseName}.log" into vcftools_log_sanity_ch
	path "${filtvcflog.baseName[0..-5]}.OK.vcf.gz" optional true into vcftools_ok_vcf_ch
	
	"""
	logstats.sh $logfile $allvcflog $filtvcflog $min_contig_length $min_filt_contig_length > ${logfile.baseName}.log
	"""
	
}
	
	
process gatkSiteFilter {

	// Perform GATK site filters on VCFs
	// Modified from RatesTools 0.5.15: Armstrong & Campana 2023
	
	publishDir "$params.outdir/04_GATKFilterChrVCFs", mode: 'copy', pattern: '*gatk.vcf.gz'
	
	input:
	path vcftools_vcf from vcftools_ok_vcf_ch
	path refseq from params.refseq
	path fai from params.refseq_fai
	path dict from params.refseq_dict
	val site_filters from params.gatk_site_filters
	
	output:
	tuple path("${vcftools_vcf.baseName[0..-17]}.gatk.tmp"), path(vcftools_vcf), path("${vcftools_vcf.baseName[0..-17]}.gatk.vcf.gz") into gatk_vcf_ch
	
	script:
	if (site_filters == "NULL")
		"""
		ln -s $vcftools_vcf ${vcftools_vcf.baseName[0..-17]}.gatk.vcf.gz
		vcftools --gzvcf $vcftools_vcf
		cp .command.log ${vcftools_vcf.baseName[0..-17]}.gatk.tmp
		"""
	else
		"""
		tabix $vcftools_vcf
		$gatk VariantFiltration -R $refseq -V $vcftools_vcf -O tmp.vcf.gz $site_filters
		$gatk SelectVariants -R $refseq -V tmp.vcf.gz -O ${vcftools_vcf.baseName[0..-17]}.gatk.vcf.gz --exclude-filtered
		rm tmp.vcf.gz
		vcftools --gzvcf ${vcftools_vcf.baseName[0..-17]}.gatk.vcf.gz
		tail .command.log > ${vcftools_vcf.baseName[0..-17]}.gatk.tmp
		"""

}

process sanityCheckLogsGatk {
	
	// Sanity check logs for GATK site filtering and remove too short contigs
	// Dummy value of 1 for min_contig_length since already evaluated and no longer accurate
	// Modified from RatesTools 0.5.15: Armstrong & Campana 2023

	input:
	tuple path(logfile), path(allvcflog), path(filtvcflog) from gatk_vcf_ch
	val min_filt_contig_length from params.min_filt_contig_length
	
	output:
	path "${logfile.baseName}.log" into gatk_sitefilt_log_sanity_ch
	path "${filtvcflog.baseName[0..-5]}.OK.vcf.gz" optional true into gatk_ok_vcf_ch
	
	"""
	logstats.sh $logfile $allvcflog $filtvcflog 1 $min_filt_contig_length > ${logfile.baseName}.log
	"""
	
}

process concatenateVCFs {

	// Concatenate VCFs using BCFtools concat
	// Attempt to sort by order in original file
	
	publishDir "$params.outdir/05_ConcatVCF", mode: 'copy'
	
	input:
	path chrfile from chrfile_ch2
	path "*" from gatk_ok_vcf_ch.collect()
	val stem from params.stem
	
	output:
	path "${stem}.all.vcf.gz" into concat_vcf_ch
	
	"""
	#!/usr/bin/env bash
	fileline=''
	readarray -t chrs < $chrfile
	len=\${#chrs[@]}
	for ((i=0; i<\$len; i++)); do
		if [ -f ${stem}_\${chrs[\$i]}.gatk.OK.vcf.gz ]; then
			fileline+="  ${stem}_\${chrs[\$i]}.gatk.OK.vcf.gz"
		fi
	done
	bcftools concat -O v -o ${stem}.all.vcf.gz\$fileline
	"""

}

process filterMappability {

	// Remove sites with mappability < 1.0 using BEDtools
	
	publishDir "$params.outdir/06_MapFiltVCF", mode: 'copy'
	
	input:
	path vcf from concat_vcf_ch
	path bed from genmap_ch
	val stem from params.stem
	
	output:
	path "${stem}.map.vcf.gz" into map_vcf_ch
	
	"""
	bedtools subtract -a $vcf -b $bed -header | gzip > ${stem}.map.vcf.gz
	"""

}

process snpRelate {

	// Estimate kinship using SNPRelate
	
	publishDir "$params.outdir/07_SNPRelate", mode: 'copy'
	
	input:
	path vcf from map_vcf_ch
	val stem from params.stem
	val snprelate_opts from params.snprelate_opts
	
	output:
	path "${stem}.gds"
	path "${stem}*.csv"
	path "${stem}.snprelate.log"
	
	"""
	#!/usr/bin/env Rscript
	library("SNPRelate")
	source(system("which kinshipUtils.R", intern = TRUE))
	snpgdsVCF2GDS(Sys.readlink(\'$vcf\'), \'${stem}.gds\', method = "biallelic.only")
	snps <- snpgdsOpen(\'${stem}.gds\')
	pruned <- snpgdsLDpruning(snps, $snprelate_opts)
	bootstrapped <- bootstrap.kinship(snps, ibdmethod = "MLE", $snprelate_opts)
	write.kinship.matrix(bootstrapped, meanfile = \"${stem}_bootstrap_meanvalues.csv\", \"${stem}_random_kinship_CI.csv\")
	system(\"cp .command.log ${stem}.snprelate.log\")
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