rule quant:
	input:
		"results_{ref}/star/{raw}Aligned.toTranscriptome.out.bam"
	output:
		"results_{ref}/salmon/{raw}/quant.sf"
	params:
		tfa = config["REF"]["TFA"]
	threads:
		64
	shell:
		"""
		salmon quant --threads {threads} \
			-t {params.tfa} \
			-l A \
			-o "results_{wildcards.ref}/salmon/{wildcards.raw}" \
			-a {input}
		"""