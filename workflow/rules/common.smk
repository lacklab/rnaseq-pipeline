import pandas as pd
#import yaml
#from pathlib import Path
#config = yaml.safe_load(Path("config/config.yaml").read_text())

# `samples` includes 'IP/CONTROL' pairs 
samples = pd.read_table(config["SAMPLES"])

# `units` includes all files that need to be preprocessed
units = pd.read_table(config["UNITS"])
units["Raw"] = units["Name"] + "_" + units["Unit"].astype(str)

ref = config["REF"]["NAME"]
# >>> utils >>>
def get_lib(wildcards):
	return units.loc[units["Name"] == wildcards.raw, "Library"].unique()[0]

def get_units(wildcards):
	return units.loc[units["Name"] == wildcards.raw, "Raw"].unique()

def get_fq1(wildcards):
	return units.loc[units["Name"] == wildcards.raw, "Fastq1"].unique()[0]

def get_fq2(wildcards):
	return units.loc[units["Name"] == wildcards.raw, "Fastq2"].unique()[0]

# <<< utils <<<



def get_fqs(wildcards):
	lib = get_lib(wildcards)
	if lib == "Single":
		return get_fq1(wildcards)
	elif lib == "Paired":
		return get_fq1(wildcards), get_fq2(wildcards)


def get_trimmed_fqs(wildcards):
	lib = get_lib(wildcards)
	if lib == "Single":
		return f"trimmed-data/{wildcards.raw}_1.fq.gz"
	elif lib == "Paired":
		return f"trimmed-data/{wildcards.raw}_1.fq.gz", f"trimmed-data/{wildcards.raw}_2.fq.gz"




outputs = []

if config["OUTPUT"]["RUN"]["QC"]:
	outputs += ["qc/multiqc_report.html"]

if config["OUTPUT"]["RUN"]["QUANT"]:
	outputs += [
		f"results_{ref}/salmon/{raw}/quant.sf"
		for raw in samples["Name"]
	]
