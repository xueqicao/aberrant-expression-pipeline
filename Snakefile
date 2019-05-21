import sys
# Add the folder path for the python parsing functions to the sys.path list
sys.path.insert(0,'../genetic_diagnosis_modified/src/python') 
from config_helper import ConfigHelper

configfile: "wbuild.yaml" 
parser = ConfigHelper(config)

# Only needed for Aberrant Expression: Outrider

outrider_all_ids, outrider_filtered = parser.getOutriderIds()
config["outrider"] = outrider_all_ids
config["outrider_filtered"] = outrider_filtered

include: ".wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables
#htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"
htmlOutputPath = "Output/html"


rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html"
    output: touch("Output/all.done")

