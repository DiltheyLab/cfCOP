#!/usr/bin/env python3

# exact copies from readstats_from_sam.py

import os, sys
import re

class SamLine:
    def __init__(self, line):
        array = line.strip().split("\t")
        self.array = array
        self.query_name = array[0]
        self.flag = int(array[1])
        self.refname = array[2]
        self.pos = int(array[3])
        self.cigar = self.array[5] 
        self.resolved_flag = self.resolve_flag() 
        self.NM = -1
        if self.array[11].startswith("NM:i:"):
            self.NM = int(self.array[11].strip("NM:i:"))
               
    def __str__(self):
        return "\t".join(self.array)

    def resolve_flag(self):
        flag=self.flag
        read_features = []
        acceptable_flags = [2048,1024,512,256,128,64,32,16,8,4,2,1]
        if flag==0:
            read_features.append(flag)
        if flag<sum(acceptable_flags):
            for af in acceptable_flags:
                if flag>=af:
                    read_features.append(af)
                    flag-=af
        return read_features
    
    
    def has_flag(self, flag): # check if samline has a given flag
        if not self.resolved_flag:
            self.resolve_flag()
        if flag in self.resolved_flag:
            return True
        else:
            return False

        
# XA tag is specific to BWA...
# (chr,pos,CIGAR,NM;)
class XA:
    def __init__(self, chr, *args):
        self.chr = chr
        self.taxid = self.chr.strip().split("|kraken:taxid|")[1]
        #self.args = args
        self.pos, self.cigar, self.nm = args
        # the XA tags in yara and bwa sam are different. 
        # BWA format:
        # Chromosome, Position(with strand), CIGAR, NM
        # Yara Format:
        # chr,begin,end,strand,NM
        
    def __str__(self):
        return "\t".join(self.chr, *self.args)


class SamLineXA(SamLine):
    def __init__(self, array):
        super().__init__(array)
        line = "\t".join(self.array)
        self.xa_tags = []
        if self.refname=="*":
            self.taxid = "0"
        else:
            self.taxid = self.refname.split(" ")[0].split("|kraken:taxid|")[1]
        
        self.NM = -1
        if self.array[11].startswith("NM:i:"):
            self.NM = int(self.array[11].strip("NM:i:"))
            
        self.xa_taxids = []
        
        split_by_XA = line.split("\tXA:Z:")
        
        if len(split_by_XA) == 2:
            xa_tags = split_by_XA[1].split("\t")[0].rstrip(";").split(";")
            for tag in xa_tags:
                try:
                    args = tag.rstrip(";").split(",")
                    xa = XA(*args)
                    self.xa_tags.append(xa)
                    self.xa_taxids.append(xa.taxid)
                except(ValueError):
                    sys.stderr.write("Warning: unsupported XA tag will be ignored. %s\n"%tag)

class CIGAR:
    def __init__(self, string):
        pattern = re.compile('([MIDNSHPX=])')
        vals = pattern.split(string)[:-1]
        self.tuples = list(([int(vals[i]), vals[i + 1]] for i in range(0, len(vals), 2)))
        self.readlen = self.summarize('MISH=X')
        #self.align_len = self.summarize('MIDN=X') # 'MIDN=X'
        self.align_len_onread = self.summarize('MI=X')
        self.insertions_to_ref = self.summarize('I')
        self.deletions_to_ref = self.summarize('D')
        
    def summarize(self, patterns):
        arr = [int(pos) for pos, match in self.tuples if match in list(patterns)]
        return sum(arr) # list('MDN=X')
    