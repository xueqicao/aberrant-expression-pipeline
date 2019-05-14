import sys
# Add the folder path for the python parsing functions to the sys.path list
sys.path.insert(0,'../genetic_diagnosis_modified/Scripts/src/python') 
import configParser
from configParser import MyConfigParser

configfile: "wbuild.yaml" 
parser = MyConfigParser(config)


# Only needed for Aberrant Expression: Outrider

outrider_all_ids, outrider_filtered = parser.geOutriderFiles("RNA_Seq")
config["outrider"] = outrider_all_ids
config["outrider_filtered"] = outrider_filtered

include: ".wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables
#htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"
htmlOutputPath = "Output/html"


rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html"
    output: touch("Output/all.done")

