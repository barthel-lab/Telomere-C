samples = [sample for sample in config["data"] if "TelomereC" in config["data"][sample]]
print(samples)
# Dict: sample:int(enrichs)
enrich = {sample: list(config["data"][sample]["TelomereC"].keys()) for sample in samples}

# References
ref_fasta="/tgen_labs/barthel/references/CHM13v2/chm13v2.0.fasta"
blacklist="data/T2T.excluderanges.noTelo.bed"

#Functions to define the input
def get_TelomereC_r1(wildcards):
  return config["data"][wildcards.sample]["TelomereC"][wildcards.enrich]["R1"]

def get_TelomereC_r2(wildcards):
  return config["data"][wildcards.sample]["TelomereC"][wildcards.enrich]["R2"]

#def get_enrichs(samples, enrichs):
#    return [enrich for sample in samples for enrich in enrichs[sample]]

#def get_enrichs(wildcards):
#  return list(config["data"][wildcards.sample]["TelomereC"].keys())