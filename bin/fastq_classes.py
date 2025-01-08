#!/usr/bin/python
import sys
import gzip


class FastqLine:
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual
        
    def __str__(self):
        return "\n".join(["@" + self.name, self.seq, "+", self.qual]) + "\n"

    
def open_compressed(infilepath):
    if infilepath.endswith('.gz'):
        f = gzip.open(infilepath, 'rt') # mode 'rt' opens file in text format, default is binary!!!
        return f
    else:
        f = open(infilepath)
        return f
    
    
def sanity_check_fastq(lines):
    invalid = 0
    if len(lines) != 4:
        invalid = 1
    if len(lines[1].strip()) != len(lines[3].strip()):
        invalid = 1
    if invalid:
        sys.stderr.write("Invalid fastq. Exiting\n") 
        for i in lines:
            sys.stderr.write(i.strip() + "\n")
        exit(1)

        
def parse_fastq(infastq):
    inf = open_compressed(infastq)
    cou = 0
    linecou = 0
    out = []
    for line in inf:
        if linecou > 0:
            out.append(line.strip())
            linecou += 1
        
        if line.startswith("@") and linecou == 0:
            linecou = 1
            out.append(line.strip()[1:])
        
        if linecou == 4:
            linecou = 0
            sanity_check_fastq(out)
            yield FastqLine(out[0], out[1], out[3])
            out = []

        
    if out:
        sanity_check_fastq(out)
        yield FastqLine(out[0], out[1], out[3])
        
        

