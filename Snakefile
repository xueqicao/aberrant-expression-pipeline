### SNAKEFILE ABERRANT EXPRESSION
import os
from config_parser import ConfigHelper

## ADD tmp/ DIR
tmpdir = config["ROOT"] + '/' + config["DATASET_NAME"] + '/tmp'
config["tmpdir"] = tmpdir
if not os.path.exists(tmpdir+'/AberrantExpression'):
    os.makedirs(tmpdir+'/AberrantExpression')
    
    
#print("In ABERRANT EXPRESSION", config)
parser = ConfigHelper(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile
htmlOutputPath = config["htmlOutputPath"]
include: os.getcwd() + "/.wBuild/wBuild.snakefile" 


## Do not delete this, needed in mergeCounts.R
config["outrider_all"], _ = parser.getOutriderIds()

rule all:
    input: rules.Index.output, htmlOutputPath + "/aberrant_expression_readme.html"
    output: touch(tmpdir + "/aberrant_expression.done")

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


   
### RULEGRAPH  
### rulegraph only works without print statements. Call <snakemake produce_graphs> for producing output

## For rule rulegraph.. copy configfile in tmp file
import oyaml
with open(tmpdir + '/config.yaml', 'w') as yaml_file:
    oyaml.dump(config, yaml_file, default_flow_style=False)

rulegraph_filename = htmlOutputPath + "/" + os.path.basename(os.getcwd()) + "_rulegraph"
dag_filename = htmlOutputPath + "/" + os.path.basename(os.getcwd()) + "_dag"

rule produce_graphs:
    input:
        expand("{graph}.{fmt}", fmt=["svg", "png"], graph=[rulegraph_filename, dag_filename])

rule create_rulegraph:
    output:
        rulegraph_filename + ".dot"
    shell:
        "snakemake --configfile " + tmpdir + "/config.yaml --rulegraph > {output}"
        
        
rule create_dag:
    output:
        dag_filename + ".dot"
    shell:
        "snakemake --configfile " + tmpdir + "/config.yaml --dag > {output}"


rule render_dot:
    input:
        "{prefix}.dot"
    output:
        "{prefix}.{fmt,(png|svg)}"
    shell:
        "dot -T{wildcards.fmt} < {input} > {output}"
