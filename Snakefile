import sys
import os

# Add the folder path for the python parsing functions to the sys.path list
sys.path.insert(0,'../genetic_diagnosis_modified/src/python') 
from config_helper import ConfigHelper

parser = ConfigHelper(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile

htmlOutputPath = config["htmlOutputPath"]

# Only needed for Aberrant Expression: Outrider

outrider_all_ids, outrider_filtered = parser.getOutriderIds()
config["outrider"] = outrider_all_ids
config["outrider_filtered"] = outrider_filtered

include: os.getcwd() + "/.wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables


rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html"
    output: touch(htmlOutputPath + "/../all.done")


# test : snakemake -n /s/project/genetic_diagnosis/processed_results/mll/v29/counts/fib_ss/total_counts.Rds
#snakemake -n /s/project/genetic_diagnosis/processed_results/mll/v29/counts/leukemia/total_counts.Rds

rule count:
    input: expand(parser.getProcDataDir() + "/{annotation}/counts/{dataset}/total_counts.Rds", annotation=config["GENE_ANNOTATION_NAMES"], dataset=parser.outrider_all)

rule filter_counts:
    input: expand(parser.getProcDataDir() + "/{annotation}/counts/{dataset}/filtered_counts.Rds", annotation=config["GENE_ANNOTATION_NAMES"], dataset=parser.outrider_all)

rule outrider:
    input: expand(parser.getProcResultsDir() + "/{annotation}/outrider/{dataset}/ods.Rds", annotation=config["GENE_ANNOTATION_NAMES"], dataset=parser.outrider_filtered)

rule counting_results:
    input: expand("Output/html/Counting/{annotation}/CountingSummary_{dataset}.html", annotation=config["GENE_ANNOTATION_NAMES"], dataset=parser.outrider_all)

rule outrider_results:
    input: expand(parser.getProcResultsDir() + "/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv", annotation=config["GENE_ANNOTATION_NAMES"], dataset=parser.outrider_all)
    
# overwriting wbuild rule output
rule rulegraph:
    shell: "snakemake --rulegraph | dot -Tsvg -Grankdir=TB > {config[htmlOutputPath]}/dep.svg"
