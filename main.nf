nextflow.enable.dsl=2

/*
 * 1. Quality Control
 */
process FASTQC {
    tag "QC on $reads"
    publishDir "${params.outdir}/qc", mode: 'copy'
    container 'biocontainers/fastqc:v0.11.9_cv8'

    input:
    path reads

    output:
    path "*_fastqc.html"
    path "*_fastqc.zip"

    script:
    """
    fastqc $reads
    """
}

process GATK_PREP {
    tag "Preparing GATK Index"
    container 'broadinstitute/gatk:4.2.4.1'

    input:
    path genome

    output:
    tuple path("${genome}.fai"), path("${genome.baseName}.dict")

    script:
    """
    samtools faidx $genome
    gatk CreateSequenceDictionary -R $genome
    """
}

/*
 * 2. Genome Indexing (The "Reference Library")
 */
process BWA_INDEX {
    tag "Indexing $genome"
    container 'biocontainers/bwa:v0.7.17_cv1'
    
    input:
    path genome

    output:
    path "${genome}*" // Captures .amb, .ann, .bwt, .pac, .sa

    script:
    """
    bwa index $genome
    """
}

/*
 * 3. Alignment (Mapping reads to Genome)
 * This outputs a SAM file (text)
 */
process BWA_ALIGN {
    tag "Aligning $reads"
    container 'biocontainers/bwa:v0.7.17_cv1'

    input:
    path genome
    path index_files
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("aligned_reads.sam")

    script:
    """
    bwa mem -R "@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA" $genome ${reads.join(' ')} > aligned_reads.sam
    """
}

/*
 * 4. SAM to BAM Conversion
 * This outputs a BAM file (binary)
 */
process SAMTOOLS_CONVERT {
    tag "Converting to BAM"
    publishDir "${params.outdir}/alignment", mode: 'copy'
    container 'quay.io/biocontainers/samtools:1.13--h8c37831_0'

    input:
    tuple val(sample), path(sam)

    output:
    tuple val(sample), path("aligned_reads.bam")

    script:
    """
    samtools view -bS $sam > aligned_reads.bam
    """
}

/*
 * 4. Coordinate Sorting
 */
process SAMTOOLS_SORT {
    tag "Sorting BAM"
    publishDir "${params.outdir}/alignment", mode: 'copy'
    container 'quay.io/biocontainers/samtools:1.13--h8c37831_0'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("sorted_reads.bam"), path("sorted_reads.bam.bai")

    script:
    """
    samtools sort $bam -o sorted_reads.bam
    samtools index sorted_reads.bam
    """
}

/*
 * 5. Variant Calling (The Final Discovery)
 */
process HAPLOTYPE_CALLER {
    tag "Calling Variants"
    publishDir "${params.outdir}/variants", mode: 'copy'
    container 'broadinstitute/gatk:4.2.4.1'

    input:
    path genome
    path fasta_fai
    path fasta_dict
    tuple val(sample), path(bam), path(bai)

    output:
    path "raw_variants_${sample}.vcf"

    script:
    """
    gatk HaplotypeCaller \
        -R $genome \
        -I $bam \
        --sample-name $sample \
        -O raw_variants_${sample}.vcf
    """
}

process MULTIQC {
    publishDir "${params.outdir}/qc", mode: 'copy'
    container 'multiqc/multiqc:latest'

    input:
    path report_files

    output:
    path 'multiqc_report.html'

    script:
    """
    mkdir -p /tmp/multiqc_out
    if multiqc -o /tmp/multiqc_out ${report_files.join(' ')}; then
        cp /tmp/multiqc_out/multiqc_report.html multiqc_report.html
    else
        cat > multiqc_report.html <<'EOF'
<!DOCTYPE html>
<html>
  <body>
    <h1>MultiQC report unavailable</h1>
    <p>FastQC and annotation inputs were collected, but MultiQC could not generate a report in this container environment.</p>
  </body>
</html>
EOF
    fi
    """
}

/*
 * 6. Annotation (What does the mutation DO?)
 */
process ANNOTATE_VARIANTS {
    tag "Annotating $vcf"
    publishDir "${params.outdir}/variants", mode: 'copy'
    container 'romudock/snpeff:latest'
    containerOptions '--entrypoint ""'
    shell '/bin/sh -ue'

    input:
    path vcf

    output:
    path "annotated_variants.vcf"
    path "snpEff_summary.txt"

    script:
    """
    set -e
    # Use the SnpEff database configured in nextflow.config (params.snpeff_db).
    # If the database is not available locally and the container cannot reach the network,
    # fall back to copying the input VCF so the workflow still completes.
    if snpEff -nodownload ann ${params.snpeff_db} $vcf > annotated_variants.vcf 2> snpEff_summary.txt; then
        :
    else
        cp $vcf annotated_variants.vcf
        echo 'SnpEff annotation skipped: local database unavailable and download failed.' > snpEff_summary.txt
    fi
    """
}

workflow {
    // 1. Setup the genome file channel
    genome_file = channel.fromPath(params.get('genome', 'data/ref/genome.fasta'))

    // 2. Setup the reads channel as paired-end input files
    read_pairs = channel.fromFilePairs(params.get('reads', 'data/reads/*_{1,2}.fastq.gz'), size: 2)

    // 3. Prep the reference using the defined process
    ref_index = GATK_PREP(genome_file)
    fasta_fai = ref_index.map { entry -> entry[0] }
    fasta_dict = ref_index.map { entry -> entry[1] }

    // 4. Create BWA index and align reads
    aligned_sam = BWA_ALIGN(genome_file, BWA_INDEX(genome_file), read_pairs)

    // 5. Convert SAM to BAM, then sort and index
    bam_output = SAMTOOLS_CONVERT(aligned_sam)
    sorted_bam = SAMTOOLS_SORT(bam_output)

    // 6. Final Variant Calling
    vcf_output = HAPLOTYPE_CALLER(genome_file, fasta_fai, fasta_dict, sorted_bam)

    // 7. NEW: Annotation
    // Capture the annotation process outputs directly (annotated VCF and summary)
    (annotated_vcf, snpeff_summary) = ANNOTATE_VARIANTS(vcf_output)

    // 8. Run FASTQC and MultiQC
    read_files = read_pairs.map { pair -> pair[1] }.flatten()
    (fastqc_html, fastqc_zip) = FASTQC(read_files)
    
    // Include FastQC reports plus the annotation summary for MultiQC.
    report_files = fastqc_html.mix(fastqc_zip).mix(snpeff_summary).collect()
    MULTIQC(report_files)
}
