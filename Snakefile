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

