import sys
# Add the folder path for the python parsing functions to the sys.path list
sys.path.append('../genetic_diagnosis_modified/Scripts/src')
import configParser

configfile: "wbuild.yaml" 
parser = configParser(config)


# set config variables
#mae
#vcfs, rnas = mae_files()
#config["vcfs"] = vcfs
#config["rnas"] = rnas
#config["mae_ids"] = list(map('-'.join, zip(vcfs, rnas)))

#outrider
#outrider_all_ids, outrider_filtered = outrider_files()
#config["outrider"] = outrider_all_ids
#config["outrider_filtered"] = outrider_filtered

include: ".wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables
#htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"
htmlOutputPath = "Output/html"


rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html"
    output: touch("Output/all.done")

