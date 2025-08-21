nextflow.enable.dsl=2

include { normalize_coverage_thresholds } from './modules/local/ops.nf'

workflow {
    take:
        FASTQ,            // reads: path
        SOI,              // sequence-of-interest FASTA: path
        REFGENOME,        // primary reference FASTA: path
        GENES_BED,        // genes bed file: path
        EXTRA_REFGENOME,  // optional extra ref: path|empty
        SAMPLE_PREFIX,    // string
        OUTDIR,           // path
        COVERAGE_THRESH   // list integers or string ("0 1 5 10")

    main:
        
        ch_prefix         = Channel.value(SAMPLE_PREFIX)
        ch_outdir         = Channel.value(OUTDIR)
        ch_cov_thresholds = normalize_coverage_thresholds(COVERAGE_THRESH)

        // 01 preprocessing ---------------------------------------------------
        porechop_out = PORECHOP( FASTQ, ch_prefix )
        fastp_out    = FASTP( porechop_out.fastq, ch_prefix )
        lengths      = FASTQ_LENGTHS( fastp_out.fastq, ch_prefix )
        nanoplot     = NANOPLOT( fastp_out.fastq, ch_prefix )

        // 02 map to Ref ------------------------------------------------------
        map_ref      = MAP_MINIMAP2( fastp_out.fastq, REFGENOME, ch_prefix, 'mapToRef' )
        cov_ref      = MOSDEPTH( map_ref.bam, ch_prefix, '02_mapToRef' )
        bw_ref       = BAMCOVERAGE( map_ref.bam, ch_prefix, '02_mapToRef' )
        paf_ref      = PAF_MINIMAP2( fastp_out.fastq, REFGENOME, ch_prefix, '02_mapToRef', 'mapToRef' )

        // 03 map to SOI ------------------------------------------------------
        map_soi      = MAP_MINIMAP2( fastp_out.fastq, SOI, ch_prefix, 'mapToSOI' )
        paf_soi      = PAF_MINIMAP2( fastp_out.fastq, SOI, ch_prefix, '03_mapToSOI', 'mapToSOI' )
        soi_mapped   = EXTRACT_MAPPED( map_soi.bam, ch_prefix )
        soi_reads    = SUBSEQ_READS( fastp_out.fastq, soi_mapped.readnames, ch_prefix, 'mapToSOI' )
        sniffles_soi = SNIFFLES( map_soi.bam, ch_prefix, '03_mapToSOI', 'mapToSOI' )

        // 04 map SOI reads to Ref -------------------------------------------
        map_soi2ref      = MAP_MINIMAP2( soi_reads.fastq, REFGENOME, ch_prefix, 'mapToSOIToRef' )
        paf_soi2ref      = PAF_MINIMAP2( soi_reads.fastq, REFGENOME, ch_prefix, '04_mapToSOIToRef', 'mapToSOIToRef' )
        sniffles_soi2ref = SNIFFLES( map_soi2ref.bam, ch_prefix, '04_mapToSOIToRef', 'mapToSOIToRef' )

        // 05 events from SOI->Ref BAM ----------------------------------------
        mos_soi2ref   = MOSDEPTH( map_soi2ref.bam, ch_prefix, '05_events', 'mapToSOIToRef' )
        events        = EVENTS_PER_COV( mos_soi2ref.per_base, GENES_BED, ch_cov_thresholds, ch_prefix )

        // reads that map to both SOI->Ref (non-empty detection happens inside)
        both_mapped      = EXTRACT_MAPPED( map_soi2ref.bam, ch_prefix, emit_reads_txt: true, tag: 'mapToSOIToRef' )
        both_reads_fq    = SUBSEQ_READS( fastp_out.fastq, both_mapped.readnames, ch_prefix, 'mapToSOIToRef' )

        // 06 optional extra reference ----------------------------------------
        extra_map     = OPTIONAL_MAP_EXTRA( fastp_out.fastq, EXTRA_REFGENOME, ch_prefix )

        // Summary delivery ----------------------------------------------------
        summary = SUMMARY(
                    ch_prefix, REFGENOME, EXTRA_REFGENOME,
                    porechop_out.summary, fastp_out.summary,
                    lengths.hist_all, lengths.hist_mid_80, lengths.hist_low_80,
                    nanoplot.stats,
                    map_ref.stats, cov_ref.summary, bw_ref.bw,
                    soi_mapped.bam_delivery, soi_reads.stats,
                    map_soi2ref.stats, sniffles_soi2ref.summary,
                    events.ev_3x, events.ev_5x,
                    both_reads_fq.stats,
                    extra_map.summary
                  )

        emit:
            summary_file = summary
}

// -------------------------------
// Processes (local DSL2 modules)
// -------------------------------

process PORECHOP {
    tag { prefix }
    publishDir "${params.outdir ?: '.'}/${prefix}/01_preprocessing", mode: 'copy'

    input:
    tuple path(fq), val(prefix)

    output:
    tuple val(prefix), path("${prefix}.trimm.fastq.gz"), path("${prefix}.Porechop.summary.txt") into porechop_ch

    script:
    """
    Pomoxis porechop \
      -i ${fq} \
      -o ${prefix}.trimm.fastq.gz \
      --threads ${task.cpus} &> ${prefix}.Porechop.summary.txt
    """
}

process FASTP {
    tag { prefix }
    cpus 4
    publishDir "${params.outdir ?: '.'}/${prefix}/01_preprocessing", mode: 'copy'

    input:
    tuple val(prefix), path(fq_in), path(porechop_summary)

    output:
    tuple val(prefix), path("${prefix}.trimm.filter.fastq.gz"),
                    path("${prefix}.fastpFiltering.summary.txt"),
                    path("${prefix}.trimm.filter.html")

    script:
    """
    fastp \
      -i ${fq_in} \
      -o ${prefix}.trimm.filter.fastq.gz \
      --qualified_quality_phred 15 \
      --unqualified_percent_limit 40 \
      --cut_front --cut_tail --cut_mean_quality 20 \
      --length_required 500 \
      -h ${prefix}.trimm.filter.html \
      -j ${prefix}.trimm.filter.json &> ${prefix}.fastpFiltering.summary.txt
    """
}

process FASTQ_LENGTHS {
    tag { prefix }
    publishDir "${params.outdir ?: '.'}/${prefix}/01_preprocessing", mode: 'copy'

    input:
    tuple val(prefix), path(fq)

    output:
    tuple val(prefix), path("${prefix}.trimm.filter.fastq.readLength"),
                    path("${prefix}.hist_all.txt"),
                    path("${prefix}.hist_mid80.txt"),
                    path("${prefix}.hist_low80.txt")

    script:
    """
    # extract read lengths
    zcat ${fq} | cut -d" " -f1 | seqkit fx2tab -nl - | cut -f2 > ${prefix}.trimm.filter.fastq.readLength
    python3 ${projectDir}/bin/fastq_lengths.py ${prefix}.trimm.filter.fastq.readLength \
      --out-all ${prefix}.hist_all.txt \
      --out-mid80 ${prefix}.hist_mid80.txt \
      --out-low80 ${prefix}.hist_low80.txt
    """
}

process NANOPLOT {
    tag { prefix }
    cpus 4
    publishDir "${params.outdir ?: '.'}/${prefix}/01_preprocessing/${prefix}.NanoPlot", mode: 'copy'

    input:
    tuple val(prefix), path(fq)

    output:
    tuple val(prefix), path('NanoStats.txt'), path('NanoPlot-report.html')

    script:
    """
    NanoTools NanoPlot -t ${task.cpus} --fastq ${fq} --loglength --N50 --title "${prefix} trimm filter" -o .
    """
}

process MAP_MINIMAP2 {
    tag { "${prefix}:${label}" }
    cpus 8
    publishDir "${params.outdir ?: '.'}/${prefix}/${folder}", mode: 'copy'

    input:
    tuple val(prefix), path(fq), path(ref), val(label)
    val folder

    output:
    tuple val(prefix), val(label), path("${prefix}.${label}.bam"), path("${prefix}.${label}.stats"), path("${prefix}.${label}.flagstat")

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont ${ref} ${fq} | \
      samtools sort -@ ${task.cpus} -o ${prefix}.${label}.bam -T ${prefix}.${label}.tmp
    samtools index ${prefix}.${label}.bam
    samtools stats ${prefix}.${label}.bam > ${prefix}.${label}.stats
    samtools flagstat ${prefix}.${label}.bam > ${prefix}.${label}.flagstat
    """
}

process MOSDEPTH {
    tag { prefix }
    cpus 4
    publishDir "${params.outdir ?: '.'}/${prefix}/${folder}", mode: 'copy'

    input:
    tuple val(prefix), val(label), path(bam)
    val folder

    output:
    tuple val(prefix), val(label), path("${prefix}.${label}.mosdepth.per-base.bed.gz"), path("${prefix}.${label}.mosdepth.summary.txt")

    script:
    """
    mosdepth ${prefix}.${label}.mosdepth --fast-mode ${bam}
    """
}

process BAMCOVERAGE {
    tag { prefix }
    cpus 4
    publishDir "${params.outdir ?: '.'}/${prefix}/${folder}", mode: 'copy'

    input:
    tuple val(prefix), val(label), path(bam)
    val folder

    output:
    tuple val(prefix), val(label), path("${prefix}.${label}.bw")

    script:
    """
    deeptools bamCoverage -p ${task.cpus} -b ${bam} -o ${prefix}.${label}.bw
    """
}

process PAF_MINIMAP2 {
    tag { prefix }
    cpus 4
    publishDir "${params.outdir ?: '.'}/${prefix}/${folder}", mode: 'copy'

    input:
    tuple val(prefix), path(fq), path(ref), val(folder), val(label)

    output:
    tuple val(prefix), val(label), path("${prefix}.${label}.PAF")

    script:
    """
    minimap2 -t ${task.cpus} -x map-ont ${ref} ${fq} > ${prefix}.${label}.PAF
    """
}

process EXTRACT_MAPPED {
    tag { prefix }
    cpus 2
    publishDir "${params.outdir ?: '.'}/${prefix}/${(params.tag ? params.tag : '03_mapToSOI')}", mode: 'copy'

    input:
    tuple val(prefix), val(label), path(bam)
    val emit_reads_txt false
    val tag 'mapToSOI'

    output:
    tuple val(prefix), path("${prefix}.mapped.${tag}.bam"), path("${prefix}.${tag}.Readnames.txt")

    script:
    """
    samtools view -bF 4 ${bam} > ${prefix}.mapped.${tag}.bam
    samtools index ${prefix}.mapped.${tag}.bam
    samtools view ${prefix}.mapped.${tag}.bam | cut -f1 | sort -u > ${prefix}.${tag}.Readnames.txt
    """
}

process SUBSEQ_READS {
    tag { prefix }
    cpus 2
    publishDir "${params.outdir ?: '.'}/${prefix}/${folder}", mode: 'copy'

    input:
    tuple val(prefix), path(fq_all), path(readnames)
    val folder

    output:
    tuple val(prefix), path("${prefix}.${folder}.Reads.fastq.gz"), path("${prefix}.${folder}.Reads.stats.txt")

    script:
    """
    # if empty, emit an empty gzip to keep the DAG alive
    bash ${projectDir}/bin/maybe_emit.sh ${readnames} ${prefix}.${folder}.Reads.fastq.gz <<'EOF'
    seqtk subseq ${fq_all} ${readnames} | gzip -c > ${prefix}.${folder}.Reads.fastq.gz
    EOF
    seqkit stats --all ${prefix}.${folder}.Reads.fastq.gz > ${prefix}.${folder}.Reads.stats.txt || true
    """
}

process SNIFFLES {
    tag { "${prefix}:${label}" }
    cpus 4
    publishDir "${params.outdir ?: '.'}/${prefix}/${folder}", mode: 'copy'

    input:
    tuple val(prefix), val(label), path(bam)
    val folder

    output:
    tuple val(prefix), val(label), path("${prefix}.${label}.vcf"), path("${prefix}.${label}.Sniffles2.txt")

    script:
    """
    sniffles \
      --input ${bam} \
      --vcf ${prefix}.${label}.vcf \
      --minsupport 1 --minsvlen 20 \
      --qc-output-all --long-ins-length 200000 --long-del-length 200000 \
      --detect-large-ins True --output-rnames &> ${prefix}.${label}.Sniffles2.txt
    """
}

process EVENTS_PER_COV {
    tag { "${prefix}:${cov}" }
    cpus 2
    publishDir "${params.outdir ?: '.'}/${prefix}/05_events", mode: 'copy'

    input:
    tuple val(prefix), val(label), path(per_base_bed_gz)
    path genes_bed
    val cov

    output:
    tuple val(prefix), val(cov), path("${prefix}.events.Cov_${cov}.bed"), path("${prefix}.events.Cov_${cov}.txt") into events_ch

    script:
    """
    zcat ${per_base_bed_gz} | awk -v COV=${cov//X/} '$4>COV {print $0}' | \
      bedtools merge -i - -d 1 -c 4,4,4 -o sum,count,mean > ${prefix}.events.Cov_${cov}.bed

    echo -e "#Chrom\tStart\tEnd\tAvgCov\tSize\tGene" > ${prefix}.events.Cov_${cov}.txt
    cat ${prefix}.events.Cov_${cov}.bed | cut -f1-3,6 | awk '{print $0"\t"$3-$2}' | \
      bedtools intersect -a - -b ${genes_bed} -wao | cut -f1-5,9 >> ${prefix}.events.Cov_${cov}.txt
    """
}

process OPTIONAL_MAP_EXTRA {
    tag { prefix }
    cpus 8
    publishDir "${params.outdir ?: '.'}/${prefix}/06_extra", mode: 'copy'

    input:
    tuple val(prefix), path(fq), val(extra_ref)

    output:
    tuple val(prefix), path('extra.map.txt') optional true

    when:
    extra_ref && extra_ref.toString() != ''

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont ${extra_ref} ${fq} | \
      samtools sort -@ ${task.cpus} -o ${prefix}.mapToExtraRef.bam -T ${prefix}.mapToExtraRef.tmp
    samtools index ${prefix}.mapToExtraRef.bam
    samtools stats ${prefix}.mapToExtraRef.bam > ${prefix}.mapToExtraRef.stats
    samtools flagstat ${prefix}.mapToExtraRef.bam > ${prefix}.mapToExtraRef.flagstat
    mosdepth ${prefix}.mapToExtraRef.mosdepth --fast-mode ${prefix}.mapToExtraRef.bam
    deeptools bamCoverage -p ${task.cpus} -b ${prefix}.mapToExtraRef.bam -o ${prefix}.mapToExtraRef.bw
    minimap2 -t ${task.cpus} -x map-ont ${extra_ref} ${fq} > ${prefix}.mapToExtraRef.PAF
    sniffles --input ${prefix}.mapToExtraRef.bam --vcf ${prefix}.mapToExtraRef.vcf \
      --minsupport 1 --minsvlen 20 --qc-output-all --long-ins-length 200000 \
      --long-del-length 200000 --detect-large-ins True --output-rnames &> ${prefix}.mapToExtraRef.Sniffles2.txt

    echo "EXTRA_REF done" > extra.map.txt
    """
}

process SUMMARY {
    tag { prefix }
    publishDir "${params.outdir ?: '.'}/${prefix}/${prefix}-Delivery", mode: 'copy'

    input:
    val prefix
    path ref
    val extra_ref
    path porechop_summary
    path fastp_summary
    path hist_all
    path hist_mid
    path hist_low
    path nanoplot_stats
    path mapref_stats
    path mosdepth_summary
    path bigwig
    path soi_bam
    path soi_reads_stats
    path soi2ref_stats
    path sniffles_soi2ref_summary
    path ev3x
    path ev5x
    path both_reads_stats
    path extra_summary optional true

    output:
    path '00_summary.txt'

    script:
    """
    {
      echo "Summary Info for sample ${prefix}"
      echo "\nPorechop Results:"; grep -E "reads had adapters trimmed|split based on middle" ${porechop_summary} || true
      echo "~.~.~.~.~"
      echo "Fastp Results:"; grep -B1 "total reads" ${fastp_summary} || true
      echo "~.~.~.~.~"
      echo "Read lengths:"; echo "All trimmed+filtered:"; cat ${hist_all}
      echo "80% excluding top and tail 10%:"; cat ${hist_mid}
      echo "80% excluding top 20% (longest):"; cat ${hist_low}
      echo "~.~.~.~.~"
      echo "Trimmed & Filter Reads Info:"; head -n 9 ${nanoplot_stats} || true
      echo "~.~.~.~.~"
      echo "Mapping to ${ref}"; echo -n "Reads Mapping To Reference: "; grep 'reads mapped:' ${mapref_stats} | cut -f3
      cat ${mosdepth_summary}
      echo "~.~.~.~.~"
      echo "Reads Mapping To SOI: "; samtools view -c ${soi_bam} || true
      echo "Reads Mapping To SOI stats:"; cat ${soi_reads_stats} || true
      echo "~.~.~.~.~"
      echo "Reads with SOI mapping to reference:"; grep 'reads mapped:' ${soi2ref_stats} | cut -f3
      echo "SVs for reads with SOI mapping to Ref:"; grep 'Wrote' ${sniffles_soi2ref_summary} || true
      echo "~.~.~.~.~"
      echo 'Events at >3X:'; cat ${ev3x} || true
      echo "~.~"
      echo 'Events at >5X:'; cat ${ev5x} || true
      echo "~.~.~.~.~"
      if [ -s ${both_reads_stats} ]; then echo "Reads Mapping To SOI AND reference stats:"; cat ${both_reads_stats}; fi
      if [ -n "${extra_ref}" ] && [ -e extra.map.txt ]; then
        echo "Mapping to EXTRAREFGENOME: ${extra_ref}"
        grep 'Wrote' ${prefix}.mapToExtraRef.Sniffles2.txt || true
      else
        echo "EXTRAREFGENOME not provided."
      fi
    } > 00_summary.txt
    """
}

// -------------------------
// Small helper operator(s)
// -------------------------

// Convert list/string thresholds to a channel of pretty tokens like ["00X","03X","05X"]
process normalize_coverage_thresholds {
    input:
    val cov_input
    output:
    path 'cov.list'
    script:
    """
    python3 - <<'PY'
import sys
val = """${cov_input}""".strip()
if not val:
    vals = [0,3,5,10]
else:
    if isinstance(val, str):
        vals = [int(x) for x in val.replace(',', ' ').split()]
    else:
        try:
            vals = list(val)
        except Exception:
            vals = [0,3,5,10]

fmt = [f"{v:02d}X" for v in vals]
open('cov.list','w').write("\n".join(fmt))
PY
    """
}

// Map Channel from file
Channel.fromPath('cov.list').splitText().map{ it.trim() }.set{ cov_tokens }

// scatter events per cov
def EVENTS_PER_COV(per_base, genes_bed, covs, prefix){
    per_base.map{ tuple(prefix, 'mapToSOIToRef', it) }
            .combine( Channel.value(genes_bed) )
            .combine( cov_tokens )
            | EVENTS_PER_COV
}