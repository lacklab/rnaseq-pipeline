import os
import urllib.request
import pandas as pd

# Load configuration and samples
# Uncomment and use if `yaml` and `Path` are available
import yaml
from pathlib import Path
references = yaml.safe_load(Path("config/references.yaml").read_text())

# `samples` includes 'IP/CONTROL' pairs
samples = pd.read_table(config["SAMPLES"])
samples["Raw"] = samples["Name"] + "_" + samples["Unit"].astype(str)
#print(samples)



# Utility functions
def get_lib(wildcards):
    return samples.loc[samples["Raw"] == wildcards.raw, "Library"].unique()[0]

def get_fastq(wildcards, fastq_col):
	#print(wildcards.raw)
	#print(samples.loc[samples["Raw"] == wildcards.raw, fastq_col])
	return samples.loc[samples["Raw"] == wildcards.raw, fastq_col].unique()[0]


def get_fqs(wildcards):
    fq1 = get_fastq(wildcards, "Fastq1")
    lib = get_lib(wildcards)
    if "SRR" in fq1:
        return (f"sra-data/{fq1}_1.fastq.gz", f"sra-data/{fq1}_2.fastq.gz") if lib == "Paired" else f"sra-data/{fq1}_1.fastq.gz"
    fq2 = get_fastq(wildcards, "Fastq2")
    return (fq1, fq2) if lib == "Paired" else fq1




outputs = []


if config["OUTPUT"]["RUN"]["QC"]:
    outputs.append("qc/multiqc_report.html")

if config["OUTPUT"]["RUN"]["QUANT"]:
	outputs += [
		f"results_{row['Genome']}/salmon/{row['Raw']}/quant.sf"
		for i, row in samples[["Raw", 'Genome']].drop_duplicates().iterrows()
	]


#print(outputs)



#import pandas as pd
#import yaml
#from pathlib import Path
#config = yaml.safe_load(Path("config/config.yaml").read_text())
#
## `samples` includes 'IP/CONTROL' pairs 
#samples = pd.read_table(config["SAMPLES"])
#
## `units` includes all files that need to be preprocessed
#units = pd.read_table(config["UNITS"])
#units["Raw"] = units["Name"] + "_" + units["Unit"].astype(str)
#
#ref = config["OUTPUT"]["REF"]
## >>> utils >>>
#def get_lib(wildcards):
#	return units.loc[units["Name"] == wildcards.raw, "Library"].unique()[0]
#
#def get_units(wildcards):
#	return units.loc[units["Name"] == wildcards.raw, "Raw"].unique()
#
#def get_fq1(wildcards):
#	return units.loc[units["Name"] == wildcards.raw, "Fastq1"].unique()[0]
#
#def get_fq2(wildcards):
#	return units.loc[units["Name"] == wildcards.raw, "Fastq2"].unique()[0]
#
## <<< utils <<<
#
#"""
#def get_fqs(wildcards):
#	name, unit = wildcards.raw.rsplit("_",1)
#	fq1 = units.loc[units["Name"] == name,"Fastq1"].unique()[0]
#	source = str(fq1).find("SRR") != -1
#	lib = get_lib(wildcards)
#	if source:
#		srr = fq1
#		if lib == "Single":
#			return f"sra-data/{srr}_1.fastq.gz"
#		elif lib == "Paired":
#			return f"sra-data/{srr}_1.fastq.gz", f"sra-data/{srr}_2.fastq.gz"
#	else:
#		fq1 = get_fq1(wildcards)
#		fq2 = get_fq2(wildcards)
#		if lib == "Single":
#			return fq1
#		elif lib == "Paired":
#			return fq1, fq2
#
#"""
#
## >>> merge >>>
#def get_batches_1(wildcards):
#	name = wildcards.raw
#	reps = units.loc[units["Name"] == name,"Fastq1"].tolist()
#	source = reps[0].find("SRR") != -1
#	if source:
#		return expand("sra-data/{srr}_1.fastq.gz", srr=reps)
#	else:
#		reps = units.loc[units["Name"] == name,"Fastq1"].tolist()
#		return reps
#	
#def get_batches_2(wildcards):
#	name = wildcards.raw
#	reps = units.loc[units["Name"] == name,"Fastq1"].tolist()
#	source = reps[0].find("SRR") != -1
#	if source:
#		return expand("sra-data/{srr}_2.fastq.gz", srr=reps)
#	else:
#		reps = units.loc[units["Name"] == name,"Fastq2"].tolist()
#		return reps
## <<< merge <<<
#
## >>> cutadapt >>>
#def get_merged_fqs(wildcards):
#	lib = get_lib(wildcards)
#	if lib == "Single":
#		return f"merged-data/{wildcards.raw}_1.fq.gz"
#	elif lib == "Paired":
#		return f"merged-data/{wildcards.raw}_1.fq.gz", f"merged-data/{wildcards.raw}_2.fq.gz"
## <<< cutadapt <<<
#
## >>> star >>>
#def get_trimmed_fqs(wildcards):
#	lib = get_lib(wildcards)
#	if lib == "Single":
#		return f"trimmed-data/{wildcards.raw}_1.fq.gz" # Name GSM etc
#	elif lib == "Paired":
#		return f"trimmed-data/{wildcards.raw}_1.fq.gz", f"trimmed-data/{wildcards.raw}_2.fq.gz"
## <<< star <<<
#
#
#
#outputs = []
#
#if config["OUTPUT"]["RUN"]["QC"]:
#	outputs += ["qc/multiqc_report.html"]
#
#if config["OUTPUT"]["RUN"]["QUANT"]:
#	outputs += [
#		f"results_{ref}/salmon/{raw}/quant.sf"
#		for raw in samples["Name"]
#	]
