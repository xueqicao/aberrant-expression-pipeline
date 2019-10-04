### SNAKEFILE ABERRANT EXPRESSION
import os
import drop

parser = drop.config(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile
include: config['wBuildPath'] + "/wBuild.snakefile"

tmpdir = os.path.join(config["ROOT"], 'tmp')
config["tmpdir"] = tmpdir
if not os.path.exists(tmpdir+'/AberrantExpression'):
    os.makedirs(tmpdir+'/AberrantExpression')
# remove dummy files if they exist
done = tmpdir + "/AE.done"
if os.path.exists(done):
    os.remove(done)

AE_ROOT = pathlib.Path(drop.__file__).parent / "modules/aberrant-expression-pipeline"

# get group subsets
config['outrider_all'] = parser.outrider_all
config['outrider_filtered'] = parser.outrider_filtered

rule all:
    input: 
        rules.Index.output, config["htmlOutputPath"] + "/aberrant_expression_readme.html",
        expand(
            parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv",
            annotation=list(config["GENE_ANNOTATION"].keys()),
            dataset=parser.outrider_filtered
        )
    output: touch(done)

rule read_count_qc:
    input:
        bam_files = lambda wildcards: parser.getFilePaths(group=wildcards.dataset, ids_by_group=config["outrider_all"], assay='RNA_ASSAY'),
        ucsc2ncbi = AE_ROOT / "resource/chr_UCSC_NCBI.txt",
        script = AE_ROOT / "Scripts/Counting/bamfile_coverage.sh"
    output:
        qc = parser.getProcDataDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/qc.tsv"
    params:
        sample_ids = lambda wildcards: parser.outrider_all[wildcards.dataset]
    shell:
        "{input.script} {input.ucsc2ncbi} {output.qc} {params.sample_ids} {input.bam_files}"

