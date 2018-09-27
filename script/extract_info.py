# This script creates the directory tree for a test run from an index file. The
# format of the index file matches the folder structure of the dataset by
# listing for each entry the APPLICATION_DIRECTORY, the MODEL_DIRECTORY and the
# MODEL_NAME separated by commas. Eg:
#
#  Photoscan , scanned-thing , thing.obj
#
# This script creates a root folder called run_TIMESTAMP and then, for each line
#  1. cds to run_TIMESTAMP/APP_DIR/MODEL_DIR
#  2. runs the executable on dataset_dir/APP_DIR/MODEL_DIR/MODEL_NAME
#  3. cds back to the invocation dir (3 levels up)

# relevant paths
index_path = ""
dataset_path = ""

import os
import os.path
import subprocess
from datetime import datetime
import json

class ParsingError(Exception): pass

def extract_fields(entry):
    fields = entry.split(',')
    if (len(fields) != 3):
        raise ParsingError(entry)
    return fields[0].strip(), fields[1].strip(), fields[2].strip()


def process_stat_file(stats_file_entry):
    with open(stats_file_entry, 'r') as stats_file:
        stats_dict = json.load(stats_file)
        row = stats_dict["mesh"]
        for i in range(32):
            texture_i = "texture_" + str(i)
            if texture_i in stats_dict:
                row += " , " + str(stats_dict[texture_i][0]["rw"]) + " , " + str(stats_dict[texture_i][0]["rh"])
            else:
                row += " , , "
        print row

# absolute path to index
if not os.path.isabs(index_path):
    index_path = os.path.abspath(index_path)

# absolute path to dataset
if not os.path.isabs(dataset_path):
    dataset_path = os.path.abspath(dataset_path)

with open(index_path, 'r') as f:
    for entry in f:

        try:
            app_dir, model_dir, model_name = extract_fields(entry)
        except ParsingError as e:
            print "Error while parsing entry at line ", current_line
            continue

        path_to_dir_stat = os.path.join(dataset_path, app_dir, model_dir);

        count = 0;
        for entry in os.listdir(path_to_dir_stat):
            if entry.endswith(".json"):
                process_stat_file(os.path.join(path_to_dir_stat, entry))
                count += 1

        if count == 0:
            print "STAT FILE MISSING AT " + path_to_dir_stat

