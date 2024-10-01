#!/usr/bin/nextflow

nextflow.enable.dsl = 2

process get_header {
    label 'process_medium'

    input:
    path vcf
	
    output:
    path 'header.vcf'

    """
	bcftools view -h --threads $task.cpus ${vcf} > header.vcf
    """
}

process get_vcf_plain_text {
    label 'process_medium'

    input:
    path vcf
	
    output:
    path 'out.vcf'

    """
	bcftools view -H --threads $task.cpus ${vcf} | grep "PASS" > out.vcf
    """
}

process split_vcf {
    label 'process_medium'

    input:
    path vcf
	
    output:
    path 'split_vcf.*'

	//split --number=${params.core} ${vcf} split_vcf.
    """
	split --lines=1000 ${vcf} split_vcf.
    """
}


process create_vcf {
    label 'process_medium'

    input:
    path vcf
	path header
	
    output:
    tuple path("*.gz") , path("*.gz.tbi")

    """
	FILENAME=\$(echo \"${vcf}.vcf\")
	cat ${header} ${vcf} > \$FILENAME
	bcftools view --threads $task.cpus --output-type "z" --threads ${task.cpus} -o \${FILENAME}.gz \$FILENAME
	bcftools index --tbi --threads $task.cpus \${FILENAME}.gz
	rm \$FILENAME
    """
}

process get_fst {
    label 'process_low'
	publishDir "${params.results_dir}/fst_csv", mode: 'copy', pattern: '*.csv'
	
    input:
    tuple path(vcf) , path(tbi)
	
    output:
    path '*.csv'

    """
	FILENAME=\$(echo \"${vcf}.csv\")
	calculate_fst.py --in_vcf ${vcf} --permut_count ${params.permut} --out_csv \$FILENAME
    """
}

// ========
// WORKFLOW
// ========


workflow {
	get_header( params.vcf )
	get_vcf_plain_text( params.vcf )
	chn_split_vcf = split_vcf( get_vcf_plain_text.out )
		| flatten
	create_vcf( chn_split_vcf , get_header.out )
	get_fst( create_vcf.out )
}

