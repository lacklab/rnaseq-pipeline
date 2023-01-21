# TODO: Trimmed outputs gives error while aligning 
# https://github.com/alexdobin/STAR/issues/1732

rule trim_fq:
	input:
		get_fqs
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
			cutadapt -a {params.FWD} -o {output.r1} {input}
			""")
		elif lib == "Paired":
			shell("""
			cutadapt -a {params.FWD} -A {params.REV} -o {output.r1} -p {output.r2} {input}
			""")





