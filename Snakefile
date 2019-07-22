### SNAKEFILE ABERRANT EXPRESSION

import os
from config_parser import ConfigHelper

parser = ConfigHelper(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile
htmlOutputPath = config["htmlOutputPath"]
include: os.getcwd() + "/.wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables
# create temporary folder
if not os.path.exists('tmp'):
    os.makedirs('tmp')

# OUTRIDER IDs
outrider_all_ids, outrider_filtered = parser.getOutriderIds()
#config["outrider"] = outrider_all_ids
#config["outrider_filtered"] = outrider_filtered
#print(config["outrider"])

rule all:
    input: rules.Index.output
    output: touch("Output/all.done")

rule count:
    input: htmlOutputPath + "/Scripts_Counting_Overview.html"

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
