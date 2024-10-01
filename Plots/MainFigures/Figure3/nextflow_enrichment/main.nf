#!/usr/bin/nextflow

nextflow.enable.dsl = 2

afbins_ch = Channel.of("ultrarare", "rare", "common")
iter_ch = Channel.of(1..params.iterations)

process splitIntoAlleleFrequencyBins {
    label 'process_low'
    publishDir "${params.intermediate_dir}/afbins-vcfs", mode: 'copy'

    input:
    path vcf
    each afbin

    output:
    tuple val(afbin), path("*.bed")

    """
    FILENAME=\$(echo \"${vcf}\" | rev | cut -f1 -d/ | rev | sed \"s/\\.vcf.*/.${afbin}.bed/g\")

    case \"${afbin}\" in
        ultrarare)
            FILTER='INFO/AF < 0.001'
            ;;
        rare)
            FILTER='INFO/AF >= 0.001 & INFO/AF < 0.01'
            ;;
        common)
            FILTER='INFO/AF >= 0.01'
            ;;
    esac

    bcftools filter --threads $task.cpus -i \"\$FILTER\" -o tmp.vcf ${vcf}
	vcf_to_bed.sh tmp.vcf > \$FILENAME
	rm tmp.vcf
    """
}

process get_bedFile_of_interest {

    input:
    path gtf
    val interest_item

    output:
    path "*.bed.gz"
	
	script:
	    if( interest_item == 'exon' )
			"""
			zcat ${gtf} > ${interest_item}.gtf
			../../../bin/gfftobed -e --GTF ${interest_item}.gtf | gzip - > ${interest_item}.bed.gz
			rm ${interest_item}.gtf
			"""
		else if ( interest_item == 'gene' )
			"""
			zcat ${gtf} > ${interest_item}.gtf
			../../../bin/gfftobed -g --GTF ${interest_item}.gtf | gzip - > ${interest_item}.bed.gz
			rm ${interest_item}.gtf
			"""
		else if ( interest_item == 'cds' )
			"""
			zcat ${gtf} > ${interest_item}.gtf
			../../../bin/gfftobed -c --GTF ${interest_item}.gtf | gzip - > ${interest_item}.bed.gz
			rm ${interest_item}.gtf
			"""
		else if ( interest_item == 'mrna' )
			"""
			zcat ${gtf} > ${interest_item}.gtf
			../../../bin/gfftobed -m --GTF ${interest_item}.gtf | gzip - > ${interest_item}.bed.gz
			rm ${interest_item}.gtf
			"""
		else if ( interest_item == 'tss' )
			"""
			zcat ${gtf} > ${interest_item}.gtf
			../../../bin/gfftobed -t --GTF ${interest_item}.gtf | gzip - > ${interest_item}.bed.gz
			rm ${interest_item}.gtf
			"""
		else if ( interest_item == 'start_codon' )
			"""
			zcat ${gtf} > ${interest_item}.gtf
			../../../bin/gfftobed -f 'start_codon' --GTF ${interest_item}.gtf | gzip - > ${interest_item}.bed.gz
			rm ${interest_item}.gtf
			"""
		else if ( interest_item == 'stop_codon' )
			"""
			zcat ${gtf} > ${interest_item}.gtf
			../../../bin/gfftobed -f 'stop_codon' --GTF ${interest_item}.gtf | gzip - > ${interest_item}.bed.gz
			rm ${interest_item}.gtf
			"""
}

process dropVariantsOverlappingExons {

    publishDir "${params.intermediate_dir}/main-vcf-afbins-exons-dropped", mode: 'copy'

    input:
    tuple val(afbin), path(bed)
    path exons_bed

    output:
    tuple val(afbin), path("*.noncoding.bed")

    """
    FILENAME=\$(echo \"${bed}\" | sed \"s/\\.bed.*/.noncoding.bed/g\")
    bedtools intersect -v -a ${bed} -b ${exons_bed} > \$FILENAME
    """
}


process intersect {

    publishDir "${params.results_dir}/main-vcf-summaries", mode: 'copy'
    
    input:
    tuple val(afbin), path(bed)
    path interest_bed
	val col_interest

    output:
    path "${afbin}_summary.json"

	"""
	bedtools intersect -wo -a ${bed} -b ${interest_bed} \
	| summarise_afbins_ORIGINAL.py --column ${col_interest} --input_file '-' \
	> ${afbin}_summary.json
	"""
}

process calculate_intron_intergenic_bias {

    publishDir "${params.results_dir}/main-vcf-summaries", mode: 'copy'
    
    input:
    tuple val(afbin), path(bed)
    path intron_bed

    output:
    path "*.json"

	"""
	bedtools intersect -wo -a ${bed} -b ${intron_bed} > ${afbin}_intron.bed
	bedtools intersect -v -a ${bed} -b ${intron_bed} > ${afbin}_intergenic.bed
	calculate_bias.py --intron_bed ${afbin}_intron.bed --intergenic_bed ${afbin}_intergenic.bed \
	> ${afbin}_bias.json
	rm ${afbin}_intron.bed
	rm ${afbin}_intergenic.bed
	"""
}

process shuffleAndDropExonsAndIntersect {

    publishDir "${params.results_dir}/iter-summaries", mode: 'copy', pattern: "*_summary.json"
	publishDir "${params.results_dir}/iter-bias", mode: 'copy', pattern: "*.bias.json"
	publishDir "${params.results_dir}/shuffled_bed/${afbin}", mode: 'copy', pattern: "*.shuffled.bed"

    input:
    tuple val(afbin), path(bed)
    each i
    path genome_file
    path bed_exclude
    path bed_include
	val col_of_interest
	
    output:
    path "*_summary.json"

	tuple val(afbin), path("*.shuffled.bed"), emit: shuffled_bed

	script:
	"""
	bedtools shuffle -i ${bed} -g ${genome_file} -noOverlapping -chrom > ${afbin}_iter_${i}.shuffled.bed 
	
	bedtools intersect -a ${afbin}_iter_${i}.shuffled.bed -b ${bed_exclude} -v \
	| bedtools intersect -a stdin -b ${bed_include} -wo \
	| summarise_afbins_ORIGINAL.py --input_file '-' --column ${col_of_interest} \
	> ${afbin}_iter_${i}_summary.json
	
	## Calculate intron vs intergenic bias ONLY applicable to cCRE
	#path "*.bias.json"

	#bedtools intersect -wo -a ${afbin}_iter_${i}.shuffled.bed -b ${params.intron_bed} > ${afbin}_intron.bed
	#bedtools intersect -v -a ${afbin}_iter_${i}.shuffled.bed -b ${params.intron_bed} > ${afbin}_intergenic.bed
	#calculate_bias.py --intron_bed ${afbin}_intron.bed --intergenic_bed ${afbin}_intergenic.bed \
	#> ${afbin}_iter_${i}.bias.json
	
	## Clean-up
	# rm ${afbin}_intron.bed
	# rm ${afbin}_intergenic.bed

	"""
}

process shuffleAndIntersect {

    publishDir "${params.results_dir}/iter-summaries", mode: 'copy', pattern: "*_summary.json"
	publishDir "${params.results_dir}/shuffled_bed/${afbin}", mode: 'copy', pattern: "*.shuffled.bed"

    input:
    tuple val(afbin), path(bed)
    each i
    path genome_file
    path bed_include
	val col_of_interest
	
    output:
    path "*_summary.json"
	tuple val(afbin), path("*.shuffled.bed"), emit: shuffled_bed

	script:
	"""
	bedtools shuffle -i ${bed} -g ${genome_file} -noOverlapping -chrom > ${afbin}_iter_${i}.shuffled.bed 
	bedtools intersect -a ${afbin}_iter_${i}.shuffled.bed -b ${bed_include} -wo \
	| summarise_afbins_ORIGINAL.py --input_file '-' --column ${col_of_interest} \
	> ${afbin}_iter_${i}_summary.json
	
	"""
}


// ========
// WORKFLOW
// ========


workflow {
	// Get bed file of interested region
	bed_of_interest_ch = Channel.empty()
	if ( params.bed ) {
		bed_of_interest_ch = params.bed
	}
	else {
		get_bedFile_of_interest( params.gencode_gtf , params.interest_region )
		bed_of_interest_ch = get_bedFile_of_interest.out
	}
	
	splitIntoAlleleFrequencyBins(params.vcf, afbins_ch)
    afbinsvcf_ch = splitIntoAlleleFrequencyBins.out
	
	// determine which column to summarise on in json 
	col_of_interest_ch = Channel.empty()
	if( params.ccre ) {
		col_of_interest_ch = 17
	}
	else {
		col_of_interest_ch = 11
	}
	
    // Main VCF
    mainVcf(afbinsvcf_ch , bed_of_interest_ch , col_of_interest_ch )

    // Randomised VCFs
    randomisedVcfs( afbinsvcf_ch , bed_of_interest_ch , col_of_interest_ch )
	
}


workflow mainVcf {
    take: 
		afbinsvcf_ch
		bed_ch
		col_of_interest_ch
		
	main:
		proc_afbinsvcf_ch = Channel.empty()
	    if( params.ccre ) {
		    dropVariantsOverlappingExons(afbinsvcf_ch, params.exons_bed)
			proc_afbinsvcf_ch = dropVariantsOverlappingExons.out  
        }
		else {
			proc_afbinsvcf_ch = afbinsvcf_ch
		}
		
		intersect( proc_afbinsvcf_ch , bed_ch , col_of_interest_ch )	
		
		//if( params.ccre ) {
        //    calculate_intron_intergenic_bias( proc_afbinsvcf_ch , params.intron_bed )
        //}
}


workflow randomisedVcfs{
    take: 
		afbinsvcf_ch
		bed_ch
		col_of_interest_ch

    main:
		if( params.ccre ) {
			shuffleAndDropExonsAndIntersect( afbinsvcf_ch, iter_ch, params.genome_file, params.exons_bed, params.bed , col_of_interest_ch )
        }
		else {
			shuffleAndIntersect( afbinsvcf_ch, iter_ch, params.genome_file, params.bed , col_of_interest_ch )
		}
		

}