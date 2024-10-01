#!/usr/bin/nextflow

nextflow.enable.dsl = 2

process mt2vcf {
   label 'process_medium'

    input:
    val mt
	
    output:
	path '*.npy'

    """
	mt2vcf.py --in_mt ${mt} 
    """
}

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
	path ethnicity_file
	
    output:
    path '*.csv'

    """
	FILENAME=\$(echo \"${vcf}.csv\")
	calculate_fst.py  --in_vcf ${vcf} \
	--ethnicity ${ethnicity_file} \
    --bootstrap ${params.permut} \
    --out_csv \$FILENAME

    """
}

// ========
// WORKFLOW
// ========


workflow {
	//mt2vcf( params.mt )
	get_header( params.vcf )
	get_vcf_plain_text( params.vcf )
	chn_split_vcf = split_vcf( get_vcf_plain_text.out )
		| flatten
	create_vcf( chn_split_vcf , get_header.out )
	get_fst( create_vcf.out , params.ethnicity  )
}

