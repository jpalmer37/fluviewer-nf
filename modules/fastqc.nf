process FASTQC {

    tag { sample_id }

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    //tuple val(sample_id), path("${sample_id}_R1.trim_p_fastqc.html"), path("${sample_id}_R2.trim_p_fastqc.html"), emit: html
    //tuple val(sample_id), path("${sample_id}_R1.trim_p_fastqc.zip"), path("${sample_id}_R2.trim_p_fastqc.zip"), emit: zip
    tuple val(sample_id), path("${sample_id}_R1.trim_p_fastqc.zip"), path("${sample_id}_R2.trim_p_fastqc.zip"), emit: zip
    tuple val(sample_id), path("${sample_id}_fastqc_provenance.yml"), emit: provenance

    script:
    """
    mkdir -p ./tmp
    printf -- "- process_name: fastqc\\n" > ${sample_id}_fastqc_provenance.yml
    printf -- "  tool_name: fastqc\\n  tool_version: \$(fastqc --version 2>&1 | sed -n '1 p')\\n" >> ${sample_id}_fastqc_provenance.yml
    fastqc \
      --threads ${task.cpus} \
      --dir ./tmp \
      ${sample_id}_R1.trim_p.fastq.gz \
      ${sample_id}_R2.trim_p.fastq.gz 
    
    """
}
 //removed from cutadapt script
 //--json=${sample_id}.cutadapt.json 
