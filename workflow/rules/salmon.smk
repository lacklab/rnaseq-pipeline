rule quant:
	input:
		"results_{ref}/star/{raw}Aligned.toTranscriptome.sorted.bam"
	output:
		"results_{ref}/salmon/{raw}/quant.sf"
	params:
		tfa = lambda wildcards: references[wildcards.ref]["TFA"]
	conda:
		"../envs/salmon.yaml"
	threads:
		16
	shell:
		"""

		salmon quant --threads {threads} \
			-t {params.tfa} \
			-l A \
			-o "results_{wildcards.ref}/salmon/{wildcards.raw}" \
			-a {input} 
		"""



	