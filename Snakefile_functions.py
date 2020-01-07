import numpy as np
import pandas as pd
from snakemake.io import expand

# TODO: depend on feature type
def getCountFiles(parser, annotation, feature_type, dataset):
    ids = parser.outrider_ids[dataset]
    file_stump = parser.getProcDataDir()
    file_stump += f"/aberrant_expression/{annotation}/{feature_type}/counts/"
    return expand(file_stump + "{sampleID}.Rds", sampleID=ids)

def getCountSettings(parser, sampleID, feature_type, setting):
    """
    params
        feature_type: type of features. If type=dassie, the features will be 
                      counted using count mode 'union'
        setting: one of ['strand', 'count_mode', 'paired_end', 'overlap']
    """
    sample = parser.sample_annotation.query("RNA_ID == @sampleID")
    # throw error if more than 1 value is returned
    if len(sample) > 1:
        raise ValueError(f"got more than 1 setting for sampleID={sampleID}, " +
                         f"feature_type={feature_type}, setting={setting}")
    
    if setting == "strand":
        return sample["STRAND"].iloc[0].lower()
        
    elif setting == "paired_end":
        return sample["PAIRED_END"].iloc[0]
        
    elif setting == "count_mode":
        if feature_type == "dassie":
            return "Union"
        else:
            return sample["COUNT_MODE"].iloc[0].capitalize()
            
    elif setting == "overlap":
        if feature_type == "dassie":
            return np.True_
        else:
            return sample["COUNT_OVERLAPS"].iloc[0]
    else:
        raise ValueError(f"'{setting}' is not a valid setting for " 
                          + "getCountSettings()")
