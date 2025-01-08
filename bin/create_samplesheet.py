#!/usr/bin/env python3

import os
import sys
import argparse
from datetime import datetime

def parser_io():
    parser = argparse.ArgumentParser(description="searches best references from given small sam and taxstat files")
    
    parser.add_argument("infiles", help="input file(s)", nargs="+")
    optNamed = parser.add_argument_group('optional arguments')
    optNamed.add_argument("-m","--mapper", help="minimap2 or bwa", required=False)
    optNamed.add_argument("-p","--platform", help="nanopore or illumina", required=False)
    
    
    return parser.parse_args()


'''
Arguments: make choices for mapper,platform

'''

def main():
    args = parser_io()
    
    infiles = []
        
    for inf in args.infiles:
        if os.path.isfile(inf):
            infiles.append(os.path.abspath(inf))
    
    if not infiles:
        sys.stderr.write('No input files found. Exiting pipeline')
        exit(1)
    
    mapper = "minimap2"
    if args.mapper:
        mapper = args.mapper
        
    platform = "nanopore"
    if args.platform:
        platform = args.platform
    
    
    time_stamp = datetime.now()
    time_string = time_stamp.strftime("%Y_%m_%d_%H_%M")
    outname = "Samplesheet_automatic_%s.txt"%time_string
    
    header = "\t".join(["id","platform","mapper","fastq1","fastq2", "sorted_sam"])
    
    with open(outname, "w") as outf:
        outf.write(header + "\n")
        for inf in infiles:
            id = os.path.basename(inf).split(".")[0]
            outarr = [id, platform, mapper, inf, "", ""]
            outf.write("\t".join(outarr) + "\n")
        


if __name__ == "__main__":
    main()



