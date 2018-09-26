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
executable_path = ""
index_path = ""
dataset_path = ""

# extra args for the exe
executable_options = [ "--regionCount=20" ]

import os
import os.path
import subprocess
from datetime import datetime

class ParsingError(Exception): pass

def extract_fields(entry):
    fields = entry.split(',')
    if (len(fields) != 3):
        raise ParsingError(entry)
    return fields[0].strip(), fields[1].strip(), fields[2].strip()



# absolute path to executable
if not os.path.isabs(executable_path):
    executable_path = os.path.abspath(executable_path)

# absolute path to index
if not os.path.isabs(index_path):
    index_path = os.path.abspath(index_path)

# absolute path to dataset
if not os.path.isabs(dataset_path):
    dataset_path = os.path.abspath(dataset_path)

# create directory for experiments
script_time = datetime.now()
root_dir_name = "run_" + script_time.strftime("%Y-%m-%d_%H-%M")
os.mkdir(root_dir_name)
os.chdir(root_dir_name)

root_dir = os.getcwd();

with open(index_path, 'r') as f:
    current_line = 0
    for entry in f:

        try:
            app_dir, model_dir, model_name = extract_fields(entry)
        except ParsingError as e:
            print "Error while parsing entry at line ", current_line
            continue

        path_to_model_dir = os.path.join(app_dir, model_dir)
        path_to_model_name = os.path.join(path_to_model_dir, model_name)
        path_to_dataset_model = os.path.abspath(os.path.join(dataset_path, path_to_model_name))
        
        print "Processing " + path_to_dataset_model
        
        os.makedirs(path_to_model_dir)
        os.chdir(path_to_model_dir)

        log_file = open("log.txt", 'w')
        call_list = [executable_path, path_to_dataset_model] + executable_options
        subprocess.call(call_list, stdout=log_file, stderr=subprocess.STDOUT)
        log_file.close()

        os.chdir(root_dir)
        ++current_line


