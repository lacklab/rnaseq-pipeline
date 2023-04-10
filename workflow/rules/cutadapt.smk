# TODO: Trimmed outputs gives error while aligning 
# https://github.com/alexdobin/STAR/issues/1732

rule trim_fq:
	input:
		get_merged_fqs
	output:
		r1 = "trimmed-data/{raw}_1.fq.gz",
		r2 = "trimmed-data/{raw}_2.fq.gz"
	params:
		FWD = config["ADAPTER_FWD"],
		REV = config["ADAPTER_REV"]
	threads:
		64
	run:
		lib = get_lib(wildcards)
		if lib == "Single":
			shell("""
			cutadapt -m 1 -a {params.FWD} -o {output.r1} {input}
			touch {output.r2}
			""")
		elif lib == "Paired":
			shell("""
			cutadapt -m 1 -a {params.FWD} -A {params.REV} -o {output.r1} -p {output.r2} {input}
			""")

rule merge_fqs_1:
	input:
		get_batches_1
	output:
		"merged-data/{raw}_1.fq.gz"
	threads:
		32
	run:
		if str(input).find(" ") != -1:
			shell("""
			zcat {input} | gzip > {output}
			""")
		else:
			shell("""
		    cp {input} {output}
			""")

rule merge_fqs_2:
	input:
		get_batches_2
	output:
		"merged-data/{raw}_2.fq.gz"
	threads:
		32
	run:
		if str(input).find(" ") != -1:
			shell("""
			zcat {input} | gzip > {output}
			""")
		else:
			shell("""
		    cp {input} {output}
			""")