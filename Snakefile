### SNAKEFILE ABERRANT EXPRESSION
import os
import drop

parser = drop.config(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile
htmlOutputPath = config["htmlOutputPath"]
include: os.getcwd() + "/.wBuild/wBuild.snakefile" 

## ADD tmp/ DIR
tmpdir = config["ROOT"] + '/' + config["DATASET_NAME"] + '/tmp'
config["tmpdir"] = tmpdir
if not os.path.exists(tmpdir+'/AberrantExpression'):
    os.makedirs(tmpdir+'/AberrantExpression')
# remove dummy files if they exist
done = tmpdir + "/AE.done"
if os.path.exists(done):
    os.remove(done)

# get group subsets
config['outrider_all'] = parser.outrider_all
config['outrider_filtered'] = parser.outrider_filtered

rule all:
    input: 
        rules.Index.output, htmlOutputPath + "/aberrant_expression_readme.html",
        expand(parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv",
	    annotation=list(config["GENE_ANNOTATION"].keys()),
	    dataset=parser.outrider_filtered)
    output: touch(done)

rule count:
    input: expand(parser.getProcDataDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/total_counts.Rds", annotation=list(config["GENE_ANNOTATION"].keys()), dataset=parser.outrider_filtered)

rule counting_results:
    input: htmlOutputPath + "/Scripts_Counting_AllDatasets.html"

rule outrider:
    input: expand(parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds", annotation=list(config["GENE_ANNOTATION"].keys()), dataset=parser.outrider_filtered)

rule outrider_results:
    input: expand(parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv", annotation=list(config["GENE_ANNOTATION"].keys()), dataset=parser.outrider_filtered)

rule read_count_qc:
    input:
        bamfiles = lambda wildcards: parser.getFilePaths(group=wildcards.dataset, ids_by_group=config["outrider_all"], assay='RNA_ASSAY'),
    output:
        qc = parser.getProcDataDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/qc.tsv"
    params:
        sample_ids = lambda wildcards: parser.outrider_all[wildcards.dataset],
        chrNames = "|".join(expand("{chr}", chr=config["chr_names"]))
    run:
        shell(f'echo "sampleID\trecord_count" > {output.qc}')
        for i in range(len(params.sample_ids)):
            sampleID = params.sample_ids[i]
            bamfile = input.bamfiles[i]
            cmd = f'samtools idxstats {bamfile} | grep -E "^({params.chrNames})" | cut -f3 | paste -sd+ - | bc'
            shell(f'count=`{cmd}`; echo "{sampleID}\t$count" >> {output.qc}')

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
