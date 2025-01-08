#!/usr/bin/env python3


import os
import sys
import argparse
import subprocess
from datetime import datetime

import re
from collections import defaultdict

def float_0_to_1(value):
    try:
        fval = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} is not a valid float.")
    
    if fval < 0.0 or fval > 1.0:
        raise argparse.ArgumentTypeError(f"{value} is not in the range 0 to 1.")
    return fval

def parse_samplesheet(infile):
    column_dict = {}
    
    with open(infile, 'r') as file:
        header = file.readline().strip().split("\t")
        column_dict = {col: [] for col in header}
        
        for line in file:
            values = line.strip().split("\t")
            for i, col in enumerate(header):
                column_dict[col].append(values[i] if i < len(values) else "")
    
    return column_dict        

def parse_nextflow_config(file_path):
    """
    Parse a Nextflow configuration file into a nested Python dictionary.
    Args:
        file_path (str): Path to the Nextflow configuration file.
    Returns:
        dict: Parsed configuration as a nested dictionary.
    """
    config = defaultdict(dict)
    stack = [config]  # Stack to handle nested blocks
    current_block = config  # Points to the current dictionary being filled

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            # Skip empty lines or comments
            if not line or line.startswith('//') or line.startswith('#'):
                continue

            # Match labeled block start (e.g., withLabel: big_mem {)
            labeled_block_start = re.match(r"withLabel:\s*([\w\-\_]+)\s*\{", line)
            if labeled_block_start:
                label_key = f"withLabel:{labeled_block_start.group(1)}"
                if label_key in current_block:
                    raise ValueError(f"Duplicate block key '{label_key}' detected.")
                new_block = defaultdict(dict)
                current_block[label_key] = new_block
                stack.append(current_block)  # Push current block to stack
                current_block = new_block  # Move into the new block
                continue

            # Match block start (e.g., executor { or $sge {)
            block_start = re.match(r"([\w\$\-]+)\s*\{", line)
            if block_start:
                block_key = block_start.group(1)
                if block_key in current_block:
                    raise ValueError(f"Duplicate block key '{block_key}' detected.")
                new_block = defaultdict(dict)
                current_block[block_key] = new_block
                stack.append(current_block)  # Push current block to stack
                current_block = new_block  # Move into the new block
                continue

            # Match block end
            if line == "}":
                if len(stack) == 0:
                    raise ValueError("Unbalanced closing brace '}' detected in the configuration file.")
                current_block = stack.pop()  # Pop the parent block
                continue

            # Match key-value pairs (e.g., queueSize = 50 or process.executor = 'local')
            key_value = re.match(r"([\w\.\-]+)\s*=\s*(.+)", line)
            if key_value:
                key, value = key_value.groups()
                # Strip quotes from string values
                if value.startswith(("'", '"')) and value.endswith(("'", '"')):
                    value = value[1:-1]
                # Convert numeric values
                elif value.isdigit():
                    value = int(value)
                elif re.match(r"^\d+\.\d+$", value):  # Float detection
                    value = float(value)
                current_block[key] = value
                continue

            # If we reach here, the line is invalid
            raise ValueError(f"Invalid line in configuration file: {line}")

    if len(stack) > 1:
        raise ValueError("Unbalanced opening brace '{' detected in the configuration file.")

    return config

def update_nextflow_config(file_path, param, new_value):
    """
    Update a parameter in a Nextflow config file while preserving indentation.

    Args:
        file_path (str): Path to the Nextflow config file.
        param (str): The parameter name to update.
        new_value (str/int/float/None): The new value for the parameter.

    Returns:
        None
    """
    # Read the content of the file
    with open(file_path, 'r') as file:
        content = file.read()

    # Determine the pattern and replacement based on the value type
    if isinstance(new_value, str):
        if new_value.startswith('/'):  # likely a path
            pattern = rf'(^\s*{param}\s*=\s*)(null|["\']?.*?["\']?)'
            replacement = rf'\1"{new_value}"'
        else:  # regular strings
            pattern = rf'(^\s*{param}\s*=\s*)(null|["\']?.*?["\']?)'
            replacement = rf'\1"{new_value}"'
    elif new_value is None:  # Handle `null`
        pattern = rf'(^\s*{param}\s*=\s*)(null|["\']?.*?["\']?)'
        replacement = rf'\1null'
    else:  # Numeric values
        pattern = rf'(^\s*{param}\s*=\s*)(null|["\']?.*?["\']?)'
        replacement = rf'\1{new_value}'

    # Replace the value using regex
    new_content, num_replacements = re.subn(pattern, replacement, content, flags=re.MULTILINE)

    if num_replacements == 0:
        print(f"Warning: Parameter '{param}' not found in {file_path}.")
    else:
        # Write the modified content back to the file
        with open(file_path, 'w') as file:
            file.write(new_content)
        print(f"Updated '{param}' to '{new_value}' in {file_path}.")
        
        
def minimap2_db_valid(config_content):
    if not "minimap2_index" in config_content["params"]:
        return False
    index_path = config_content["params"]["minimap2_index"]
    if os.path.isfile(index_path) and index_path.endswith('.mmi'):
        return True
    return False

def bwa_db_valid(config_content):
    if "bwa_index" not in config_content["params"]:
        return False
    index_path = config_content["params"]["bwa_index"]
    extensions = [".amb",".ann",".bwt",".pac",".sa"]
    for ext in extensions:
        if not os.path.isfile(index_path + ext):
            return False
    return True

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="cfCOP")
    parser.add_argument(
        "--samplesheet", 
        help="Path to the sample sheet (CSV or TXT format)."
    )
    parser.add_argument(
        "--inputs", 
        nargs="+", 
        help="List of input files (for multiple files)."
    )
    
    parser.add_argument(
        "--mapper", 
        choices=['minimap2', 'bwa'], 
        required=False,
        help='Specify the aligner to use. Choices are: minimap2 or bwa.'
    )
    
    parser.add_argument(
        "--platform", 
        choices=['nanopore', 'illumina'], 
        required=False,
        help='Specify the sequencing platform of sample(s). Choices are: nanopore or illumina'
    )
    
    parser.add_argument(
        '--kmer_cutoff',
        type=float_0_to_1,
        required=False,
        help="Report generation: Kmer index cutoff should be a float number between 0 and 1"
    )   
    
    parser.add_argument(
        '--count_cutoff',
        type=int,
        required=False,
        help="Report generation: number of reads assigned to a taxon in order to be reported.  Should be an integer"
    )   
    
    parser.add_argument(
        '--fraction_cutoff',
        type=float_0_to_1,
        required=False,
        help="Report generation: fraction of reads (of total reads) assigned to a taxon in order to be reported. Should be a float number between 0 and 1"
    )   
    
    parser.add_argument(
        '--connect_db',
        required=False,
        help="Connecting a database. Filling the paths of the database files into the nextflow.config"
    )   
    
    
    parser.add_argument(
        '--resume',
        action='store_true',
        help="Resume the pipeline with changed parameters or for debuggung purposes. Stored unchanged processes will be cached, only changes will be re-run"
    )

    
    args = parser.parse_args()
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Validating pipeline
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # directory with pipeline main and wrapper scripts
    base_dir = os.path.dirname(__file__)
    
    main_workflow = os.path.join(base_dir, "main.nf")
    if not os.path.isfile(main_workflow):
        sys.stderr.write("Main pipeline is missing in the directory where this script is stored. No pipeline - no analysis. Bye!\n")
        sys.exit(1)
    
    base_config = os.path.join(base_dir, "nextflow.config")
    if not os.path.isfile(base_config):
        sys.stderr.write("No base configuration file found. No specifications - no analysis. Bye!\n")
        sys.exit(1)
    
    config_content = parse_nextflow_config(base_config)    
    if not "params" in config_content:
        sys.stderr.write("Base configuration file nextflow.config is missing all parameters defined in scope 'params'. Please add parameters to nextflow.config. Bye!\n")
        sys.exit(1)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Connecting database
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    connect_failed = 0
    if args.connect_db:
        connect_db = os.path.abspath(args.connect_db)
        if not os.path.isdir(args.connect_db):
            connect_failed = 1
        db_files = ["combined_db.fasta","combined_db.mmi","Filter/complexity.bed","Filter/homology_family.bed","Filter/plasmid.bed", "taxonomy"] 
        db_params = ["database_fasta","minimap2_index","complexity_filter","homology_filter","plasmid_filter","taxonomy"]
        db_paths = [os.path.join(connect_db, x) for x in db_files]
        for p in db_paths:
            if not os.path.exists(p):
                connect_failed = 1
        
        for param, value in zip(db_params, db_paths):
            update_nextflow_config(base_config, param, value)
        if connect_failed:
            sys.stderr.write("Failed to connect to the cfCOP database. Please read documentation at https://github.com/DiltheyLab/cfCOP\n")
            sys.exit(1)
        sys.stderr.write("cfCOP database connected successfully. Try running cfCOP\n")
        sys.exit(0)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Validating input
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if args.samplesheet and args.inputs:
        sys.stderr.write("You supplied both samplesheet and input samples. Please start the run with EITHER samplesheet OR input samples. Bye!\n")
        sys.exit(1)
    
    if not args.samplesheet and not args.inputs:
        print("Error: You must provide either --samplesheet or --inputs.", file=sys.stderr)
        sys.exit(1)

    if args.samplesheet and not os.path.exists(args.samplesheet):
        print(f"Error: The specified sample sheet '{args.samplesheet}' does not exist.", file=sys.stderr)
        sys.exit(1)

    if args.inputs:
        for input_file in args.inputs:
            if not os.path.exists(input_file):
                print(f"Error: The input file '{input_file}' does not exist.", file=sys.stderr)
                sys.exit(1)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Validating mappers and databases
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    db_params = ["database_fasta","complexity_filter","homology_filter","plasmid_filter","taxonomy"]
    validation_fail = 0
    for param in db_params:
        if not param in config_content["params"]:
            validation_fail = 1
            break
        if not os.path.exists(config_content["params"][param]):
            validation_fail = 1
            break                
    if validation_fail:            
        sys.stderr.write("cfCOP database is invalid. To connect cfCOP database: run_cfCOP.py --connect_db /your/path/to/cfCOP_database \nBye!\n")
        exit(1)

    acceptable_mappers = set(["minimap2", "bwa"])
    
    if args.samplesheet:
        samplesheet = parse_samplesheet(args.samplesheet)
        used_mappers = set(samplesheet["mapper"])   
        extra = used_mappers - acceptable_mappers
        if extra:
            sys.stderr.write("Unsupported mappers: " + str(extra) + ". Please check your nextflow.config. Acceptable mappers are: minimap2, bwa.\nBye!\n")
            sys.exit(1)
            
    if args.inputs:
        if args.mapper:
            used_mappers = set([args.mapper])
        else:
            if not "mapper" in config_content["params"]:
                sys.stderr.write("Please define default mapper (mapper='minimap2') in 'params' scope of base configuration file nextflow.config. Acceptable mappers are: minimap2, bwa.\nBye!\n")
                sys.exit(1)
            used_mappers = set([config_content["params"]["mapper"]])
            #print(config_content["params"]["minimap2_index"])
            extra = used_mappers - acceptable_mappers
            if extra:
                sys.stderr.write("Unsupported mappers: " + str(extra) + ". Please check your nextflow.config. Acceptable mappers are: minimap2, bwa.\nBye!\n")
                sys.exit(1)
            
    if "minimap2" in used_mappers:
        if not minimap2_db_valid(config_content):
            sys.stderr.write("Please connect minimap2 database. Write minimap2_index=/path/to/db.mmi in your nextflow.config\nBye!\n")
            exit(1)
    if "bwa" in used_mappers:
        #print(config_content["params"]["bwa_index"])
        if not bwa_db_valid(config_content):
            sys.stderr.write("Please connect bwa database. Write bwa_index=/path/to/db in your nextflow.config\nBye!\n")
            exit(1)
            
    #print("Bye!")
    #exit(0)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Constructing Nextflow command
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nextflow_command = ["nextflow", "run", main_workflow]
        
    if args.samplesheet:
        nextflow_command.extend(["--samplesheet", args.samplesheet])
    # ToDo: do not accept and revise the second argument if samplesheet is given!!
    # ToDo: implement platform parameter here. Only if input and not samplesheet
    else:
        if args.inputs:
            # Pass multiple input files as a comma-separated string
            input_files = ",".join(args.inputs)
            nextflow_command.extend(["--input_files", input_files])
        if args.mapper:
            nextflow_command.extend(["--mapper", args.mapper])
        if args.platform:
            nextflow_command.extend(["--platform", args.platform])
    if args.kmer_cutoff:
        nextflow_command.extend(["--kmer_cutoff", str(args.kmer_cutoff)])
    if args.count_cutoff:
        nextflow_command.extend(["--count_cutoff", str(args.count_cutoff)])
    if args.fraction_cutoff:
        nextflow_command.extend(["--fraction_cutoff", str(args.fraction_cutoff)])
    if args.resume:
        nextflow_command.extend(["-resume"])
    
    # Print and run the command
    print(f"Running command: {' '.join(nextflow_command)}")
    result = subprocess.run(nextflow_command)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()



