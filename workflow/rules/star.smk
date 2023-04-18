rule align:
	input:
		#get_fqs
		get_trimmed_fqs
	output:
		"results_{ref}/star/{raw}Aligned.sortedByCoord.out.bam",
		"results_{ref}/star/{raw}Aligned.toTranscriptome.out.bam"
	params:
		idx = config[f"REF_{ref}"]["STAR_IDX"],
		gtf = config[f"REF_{ref}"]["GTF"]
	threads:
		64
	shell:
		"""
		STAR --genomeDir {params.idx} \
			--readFilesIn {input} \
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
			--quantTranscriptomeBan IndelSoftclipSingleend \
			--limitBAMsortRAM 32000000000
		"""

# https://ycl6.gitbook.io/guide-to-rna-seq-analysis/raw-read-processing/mapping/alignment-based-method