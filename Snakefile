### SNAKEFILE ABERRANT EXPRESSION
import os

# Add the folder path for the python parsing functions to the sys.path list
from config_parser import ConfigHelper

parser = ConfigHelper(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile

htmlOutputPath = config["htmlOutputPath"] if (config["htmlOutputPath"] != None) else "Output/html"

# Only needed for Aberrant Expression: Outrider
#outrider_all_ids, outrider_filtered = parser.getOutriderIds()
#config["outrider"] = outrider_all_ids
#config["outrider_filtered"] = outrider_filtered

######## test
#print(outrider_filtered)
#print(parser.outrider_filtered)
#print(parser.getOutriderIds())
#import pandas as pd
#print(parser.getSampleIDs("RNA_ID"))
#print(parser.sample_file_mapping)


include: os.getcwd() + "/.wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables


rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html"
    output: touch(htmlOutputPath + "/../all.done")


rule count:
    input: expand(parser.getProcDataDir() + "/{annotation}/counts/{dataset}/total_counts.Rds", annotation=list(config["GENE_ANNOTATION"].keys()), dataset=parser.outrider_all)

rule filter_counts:
    input: expand(parser.getProcDataDir() + "/{annotation}/counts/{dataset}/filtered_counts.Rds", annotation=list(config["GENE_ANNOTATION"].keys()), dataset=parser.outrider_all)

rule outrider:
    input: expand(parser.getProcResultsDir() + "/{annotation}/outrider/{dataset}/ods.Rds", annotation=list(config["GENE_ANNOTATION"].keys()), dataset=parser.outrider_filtered)

rule counting_results:
    input: expand(htmlOutputPath + "/Counting/{annotation}/CountingSummary_{dataset}.html", annotation=list(config["GENE_ANNOTATION"].keys()), dataset=parser.outrider_all)

rule outrider_results:
    input: expand(parser.getProcResultsDir() + "/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv", annotation=list(config["GENE_ANNOTATION"].keys()), dataset=parser.outrider_all)


# overwriting wbuild rule output
rule rulegraph:
    shell: "snakemake --rulegraph | dot -Tsvg -Grankdir=TB > {config[htmlOutputPath]}/dep.svg"
