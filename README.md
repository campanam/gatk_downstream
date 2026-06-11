# gatk_downstream  
<img align="right" src="NZP-20180628-509SB_thumb.jpg">  

Michael G. Campana, 2023-2026  
Smithsonian's National Zoo & Conservation Biology Institute  

Nextflow pipeline [1] for joint-genotyping and variant filtration using the Genome Analysis Toolkit [2] and VCFtools [3]. The pipeline also removes variants from difficult-to-align regions using GenMap [4], scripts from RatesTools [5], and BEDTools [6]. Kinship estimation is performed using SNPRelate [7] (with bootstrapping from kinshipUtils [8]) and ngsRelateV2 [9].  

## Citation  
Please cite:  
Gregory, J.T., Maldonado, J.E., Brown, J.L., McInerney, N.R., Rogers, R., Campana, M.G., Prado, N.A. In prep. Integrating genomics and endocrinology to investigate hyperprolactinemia in female African savannah elephants housed in North American zoos.  

## License  
This software is licensed under the Smithsonian Institution [terms of use](https://www.si.edu/termsofuse).  

## Installation  
After installing [Nextflow](https://www.nextflow.io/) and [Conda](https://docs.conda.io/en/latest/), download the pipeline using:  
`nextflow pull campanam/gatk_downstream`  

## Configuring the Pipeline  
The `nextflow.config` file included with this repository contains a standard profile for running the pipeline locally. See the Nextflow documentation for assistance in generating a configuration profile for your computing system. The parameters you will need to provide to execute the pipeline are listed in the `params` block. These are:  
`outdir`: Path to the output directory  
`stem`: Stem for output file names  
`refseq`: Path to the genome reference sequence  
`gvcfs`: Path to directory of per-individual gVCFs (e.g as output via the [campanam/Elephants](https://github.com/campanam/Elephants) pipeline)  
`chrlist`: Path to list of chr to genotype. Set to "NULL" to use all chromosomes.  
`java_options`: String of options for Java executables.  
`vcftools_site_filters`: Site filters to pass to VCFtools. Set to "NULL" to ignore this filter.  
`gatk_site_filters`: Site filters to pass to GATK. Set to 'NULL' to ignore this filter.  
`min_contig_vars`: Minimum number of variants on a contig (before filtering) to retain in analysis  
`min_filt_contig_vars`: Minimum number of variants on a contig (after filtering) to retain in analysis  
`gm_tmpdir`: Scratch directory for GenMap indexing  
`gm_opts`: Specifies mapping parameters (as a string) for GenMap. Use cpus in the process configuration to set number of concurrent threads.  
`snprelate`: Run SNPRelate (true or false)  
`snprelate_ld`: LD-threshold for SNPRelate  
`snprelate_opts`: String of other SNPRelate options (excluding LD-threshold)  
`ngsrelate`: Run ngsRelateV2 (true or false)  
`ngsrelate_opts`: String of options for ngsRelateV2 (other than number of threads)  
`email`: Email to send completion status to. Set to "NULL" for no email.  

## Executing the Pipeline  

Execute the pipeline using the following command:  
`nextflow run campanam/gatk_downstream -r main -c <config_file.config> -profile standard`  

## References  
1. Di Tommaso, P., Chatzou, M., Floden, E.W., Prieto Barja, P., Palumbo, E., Notredame, C. (2017) Nextflow enables reproducible computational workflows. *Nat Biotechnol*, __35__, 316–319. DOI: [10.1038/nbt.3820](https://www.nature.com/articles/nbt.3820).  
2. McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., DePristo, M.A. (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res*, __20__, 1297-1303. DOI: [10.1101/gr.107524.110](https://genome.cshlp.org/content/20/9/1297.abstract).  
3. Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo, M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T., McVean, G., Durbin, R. (2011) The variant call format and VCFtools. *Bioinformatics*, __27__, 2156–2158. DOI: [10.1093/bioinformatics/btr330](https://academic.oup.com/bioinformatics/article/27/15/2156/402296).  
4. Pockrandt, C., Alzamel, M., Iliopoulos, C.S., Reinert, K. (2020) GenMap: ultra-fast computation of genome mappability. *Bioinformatics*, __36__, 3687–3692. DOI: [10.1093/bioinformatics/btaa222](https://academic.oup.com/bioinformatics/article/36/12/3687/5815974).
5. Armstrong, E.E., Campana, M.G. 2023. RatesTools: a Nextflow pipeline for detecting *de novo* germline mutations in pedigree sequence data. *Bioinformatics*, __39__, btac784. DOI: [10.1093/bioinformatics/btac784](https://doi.org/10.1093/bioinformatics/btac784).  
6. Quinlan, A.R., Hall, I.M. (2010) BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*, __26__, 841-842, DOI: [10.1093/bioinformatics/btq0333](https://academic.oup.com/bioinformatics/article/26/6/841/244688).  
7. Zheng, X., Levine, D., Shen, J., Gogarten, S.M., Laurie, C., Weir, B.S. (2012) A high-performance computing toolset for relatedness and principal component analysis of SNP data. *Bioinformatics*, __28__, 3326-3328. DOI: [10.1093/bioinformatics/bts606](https://doi.org/10.1093/bioinformatics/bts606).  
8. Cortes-Rodriguez, N., Campana, M.G., Berry, L., Faegre, S., Derrickson, S.R., Ha, R.R., Dikow, R.B., Rutz, C., Fleischer, R.C. (2019) Population genomics and structure of the critically endangered Mariana crow (*Corvus kubaryi*). *Genes*, __10__, 187. DOI: [10.3390/genes10030187](https://doi.org/10.3390/genes10030187).  
9. Hanghøj, K., Moltke, I., Andersen, P.A., Manica, A., Korneliussen, T.S. (2019) Fast and accurate relatedness estimation from high-throughput sequencing data in the presence of inbreeding. *GigaScience*, __8__, giz034. DOI: [10.1093/gigascience/giz034](https://doi.org/10.1093/gigascience/giz034).  

Image Credit: Skip Brown. 2018. Smithsonian's National Zoo & Conservation Biology Institute. Smithsonian Institution. https://www.si.edu/object/asian-elephant:nzp_NZP-20180628-509SB.  
