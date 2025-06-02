rule star_index:
	input:
		genome_fasta=lambda wildcards: references[wildcards.ref]["FA"],
		annotation_gtf=lambda wildcards: references[wildcards.ref]["GTF"],
		transcript_fasta=lambda wildcards: references[wildcards.ref]["TFA"]
	output:
		directory("results_{ref}/star_index/")
	params:
		sjdbOverhang=99  # Ideally set to read length - 1
	conda:
		"../envs/star.yaml"
	resources:
		mem_mb=700000,
		partition='long,big-mem',
		runtime=144000
	threads: 
		32
	shell:
		"""
		STAR --limitGenomeGenerateRAM 500000000000\
			--runMode genomeGenerate \
			--runThreadN {threads} \
			--genomeDir results_{wildcards.ref}/star_index/ \
			--genomeFastaFiles {input.genome_fasta} {input.transcript_fasta} \
			--sjdbGTFfile {input.annotation_gtf} \
			--sjdbOverhang {params.sjdbOverhang}
		"""

rule star:
	input:
		trimmed_fq1="trimmed/{raw}_1.trimmed.fastq.gz",
		trimmed_fq2="trimmed/{raw}_2.trimmed.fastq.gz",
		prereq = "results_{ref}/star_index/"
	output:
		"results_{ref}/star/{raw}Aligned.sortedByCoord.out.bam",
		"results_{ref}/star/{raw}Aligned.toTranscriptome.out.bam"
	params:
		idx = lambda wildcards: f"results_{wildcards.ref}/star_index/",
		gtf = lambda wildcards: references[wildcards.ref]["GTF"]
	conda:
		"../envs/star.yaml"
	threads:
		16
	shell:
		"""
		STAR --genomeDir {params.idx} \
			--readFilesIn {input.trimmed_fq1} {input.trimmed_fq2} \
			--outFileNamePrefix results_{wildcards.ref}/star/{wildcards.raw} \
			--readFilesCommand zcat \
			--runThreadN {threads} \
			--genomeLoad NoSharedMemory \
			--twopassMode Basic \
			--sjdbGTFfile {params.gtf} \
			--sjdbScore 2 \
			--sjdbOverhang 99 \
			--limitSjdbInsertNsj 1000000 \
			--outFilterMultimapNmax 20 \
			--alignSJoverhangMin 8 \
			--alignSJDBoverhangMin 1 \
			--outFilterMismatchNmax 999 \
			--outFilterMismatchNoverReadLmax 0.04 \
			--alignIntronMin 20 \
			--alignIntronMax 1000000 \
			--alignMatesGapMax 1000000 \
			--outSAMunmapped Within \
			--outFilterType BySJout \
			--outSAMattributes NH HI AS NM MD \
			--outSAMtype BAM SortedByCoordinate \
			--quantMode TranscriptomeSAM GeneCounts \
			--quantTranscriptomeSAMoutput BanSingleEnd_BanIndels_ExtendSoftclip \
			--limitBAMsortRAM 300000000000
		"""

rule star_sort:
	input:
		"results_{ref}/star/{raw}Aligned.toTranscriptome.out.bam"
	output:
		"results_{ref}/star/{raw}Aligned.toTranscriptome.sorted.bam"
	conda:
		"../envs/samtools.yaml"
	threads:
		16
	shell:
		"""
		samtools sort -@ {threads} -o {output} {input}

		samtools index {output}
		"""

# https://ycl6.gitbook.io/guide-to-rna-seq-analysis/raw-read-processing/mapping/alignment-based-method