### SNAKEFILE ABERRANT EXPRESSION
import os
import drop
import pathlib

parser = drop.config(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile
include: config['wBuildPath'] + "/wBuild.snakefile"

METHOD = 'AE'
SCRIPT_ROOT = drop.getMethodPath(METHOD, link_type='workdir')
TMP_DIR = config['tmpdir']

# get group subsets
config['outrider_all'] = parser.outrider_all
config['outrider_filtered'] = parser.outrider_filtered

rule all:
    input: 
        rules.Index.output, config["htmlOutputPath"] + "/aberrant_expression_readme.html",
        expand(
            config["htmlOutputPath"] + "/AberrantExpression/Counting/{annotation}/Summary_{dataset}.html",
            annotation=list(config["geneAnnotation"].keys()),
            dataset=parser.outrider_filtered
        ),
        expand(
            parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv",
            annotation=list(config["geneAnnotation"].keys()),
            dataset=parser.outrider_filtered
        )
    output: touch(drop.getMethodPath(METHOD, link_type='final_file', tmp_dir=TMP_DIR))

rule read_count_qc:
    input:
        bam_files = lambda wildcards: parser.getFilePaths(group=wildcards.dataset, ids_by_group=config["outrider_all"], file_type='RNA_BAM_FILE'),
        ucsc2ncbi = os.path.join(SCRIPT_ROOT, "resource", "chr_UCSC_NCBI.txt"),
        script = os.path.join(SCRIPT_ROOT, "Scripts", "Counting", "bamfile_coverage.sh")
    output:
        qc = parser.getProcDataDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/bam_coverage.tsv"
    params:
        sample_ids = lambda wildcards: parser.outrider_all[wildcards.dataset]
    shell:
        "{input.script} {input.ucsc2ncbi} {output.qc} {params.sample_ids} {input.bam_files}"


### RULEGRAPH
import oyaml

config_file = drop.getMethodPath(METHOD, link_type='config_file', tmp_dir=TMP_DIR)
rulegraph_filename = f'{config["htmlOutputPath"]}/{METHOD}_rulegraph'

with open(config_file, 'w') as yaml_file:
    oyaml.dump(config, yaml_file, default_flow_style=False)

rule produce_rulegraph:
    input:
        expand(rulegraph_filename + ".{fmt}", fmt=["svg", "png"])

rule create_graph:
    output:
        rulegraph_filename + ".dot"
    shell:
        "snakemake --configfile {config_file} --rulegraph > {output}"

rule render_dot:
    input:
        "{prefix}.dot"
    output:
        "{prefix}.{fmt,(png|svg)}"
    shell:
        "dot -T{wildcards.fmt} < {input} > {output}"
