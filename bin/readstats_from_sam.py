#!/usr/bin/env python3

import os
import sys
import argparse
from datetime import datetime
import re

# sys.path.append('/home/tyshaiev/my_scripts/refseq_db_download_and_indexing')
import taxonomy_NCBI as taxo


def parser_io():
    parser = argparse.ArgumentParser(description="this script makes (partial) LCA based on the BWA XA tags in a given sam")
    
    parser.add_argument("insam", help="input sam file")
    parser.add_argument("out", help="Path to the output file")
    parser.add_argument("taxonomy", help="NCBI taxonomy folder (names.dmp, nodes.dmp, merged.dmp neccessary)")
    parser.add_argument("--complexity", help="read names newline-delimited which did not survive complexity filtering (sDUST)")
    parser.add_argument("--homology", help="read names newline-delimited which did not survive homology filtering (self-mapping strategy)")
    parser.add_argument("--plasmid", help="read names newline-delimited which map to plasmids")
    return parser.parse_args()
    
    
# This script:
# 1. read line. state machine, per read
# 2. read xa tags, read supplement alignments
# 3. append to read stats, krakenlike format
# 4. make read taxid as LCA of all supplement and XA alignmets

'''
new to v0.2 of pipeline:
last columns are 1 or 0 for read crosses some masking in database
therefore readfile gets a header!!!

new to v0.5:
input: filter are optional arguments now
cigar string alignment length is counted on read only (ci.align_len_onread), so that alignment length is never over 100%
self.align_len_onread = self.summarize('MI=X')
'''


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


def get_to_the_root(node_dict, node, way):
    way.append(node.taxid)
    if node.taxid=="1":
        way.reverse()
        return way
    return get_to_the_root(node_dict, node_dict[node.parent], way)


def LCA(node_dict, *taxids):
    if len(set(taxids)) == 1: # because many alternative human locations
        return taxids[0]
    id_list = [get_to_the_root(node_dict, node_dict[id], []) for id in taxids]
    minlen = min([len(i) for i in id_list])
    pos = 0
    for pos in range(minlen):
        id_set = set([to_root[pos] for to_root in id_list])
        if len(id_set) > 1:
            return(id_list[0][pos-1])
    return id_list[0][pos]



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
    
    
class Stat:
    def __init__(self, map_type, ref, nm, align_len, insert, delet):
        self.map_type = map_type # M for main, X for XA-tag, S for supplementary alignment, U for unmapped
        self.nm = int(nm)
        self.align_len = int(align_len)
        self.insert = int(insert)
        self.delet = int(delet)
        if self.map_type == "U": # unmapped
            self.taxid = "0"
        else:
            self.taxid = ref.strip().split("|kraken:taxid|")[1]
    
    def __str__(self):
        #if self.map_type == "U":
        #    return "-"
        strings = [str(x) for x in [self.map_type, self.taxid, self.nm, self.align_len, self.insert, self.delet]]
        return ":".join(strings)

class StatLine:
    def __init__(self, readname):
        self.flag = "U"
        self.name = readname
        self.taxid = "0"
        self.length = [0,0] # length for forwards and  reverse read
        self.stats = [[],[]] # stats for forwards and  reverse read
    
    def __str__(self):
        if not self.length[1]: # single end
            length = str(self.length[0])
            stats = " ".join([x.__str__() for x in self.stats[0]])
            
        elif self.length[0] and self.length[1]:
            length = "|".join([str(x) for x in self.length])
            stats = []
            for st in self.stats:
                stats.append(" ".join([x.__str__() for x in st]))
            stats = " |:| ".join(stats)
            
        return "\t".join([self.flag, self.name, self.taxid, length, stats])
    
    def add_samline(self, samline):
        if samline.has_flag(256) or samline.has_flag(2048): # not primary or supplementary alignment
            map_type = "S"
        else:
            map_type = "M"
        
        if samline.has_flag(4): # read unmapped
            readlen = len(samline.array[9]) # read sequence
            stat = Stat("U", 0,0,0,0,0) 
        else:
            ci = CIGAR(samline.cigar)
            readlen = ci.readlen
            stat = Stat(map_type, samline.refname, samline.NM, ci.align_len_onread, ci.insertions_to_ref, ci.deletions_to_ref)
            self.flag = "C"
            
        where_to_put = 0
        if samline.has_flag(128): # second in pair
            where_to_put = 1
        
        
        # HERE WE SUMMARIZE LENGTH OF ALL MAIN MATCHES!! WORKAROUND FOR IMPAIRED READS, LIKE BWA-A
        if map_type == "M":
            self.length[where_to_put] += readlen
        #if not self.length[where_to_put]:
        #    self.length[where_to_put] = readlen
        #elif self.length[where_to_put] != readlen:
        #    print("something is wrong!! Read length was %d and changed to %d"%(self.length[where_to_put],readlen))
        
        self.stats[where_to_put].append(stat)
        
        for xa in samline.xa_tags:
            ci = CIGAR(xa.cigar)
            stat = Stat("X", xa.chr, xa.nm, ci.align_len_onread, ci.insertions_to_ref, ci.deletions_to_ref)
            self.stats[where_to_put].append(stat)
            
    def do_lca(self, nodedict):
        taxids = []
        for rst in self.stats:
            for st in rst:
                # currently all alignments will be considered for the lca!!!
                #if st.map_type != "S": # do not consider supplementary alignments for LCA (only X and M are considered. U if both unmapped)
                taxids.append(st.taxid)
        
        taxids = list(set(taxids))
        if len(taxids) == 1:
            self.taxid = taxids[0]
        else:
            if "0" in taxids:
                taxids.remove("0") # if one read is unmapped, only the mapped one will be considered
                if len(taxids) == 1:
                    self.taxid = taxids[0]
                    return               
            self.taxid = LCA(nodedict, *taxids)
            
class Filter:
    def __init__(self):
        self.maindict = {"complexity": "unset", "homology": "unset", "plasmid": "unset"}
        
    def add_filter(self, filepath, filtername):
        if os.path.exists(filepath):
            filter = set()
            with open(filepath) as ff:
                for line in ff:
                    readname = line.strip()
                    if readname:
                        filter.add(readname)
            
            self.maindict[filtername] = filter
                
    def filter_readname(self, readname):
        outline = []
        for filt in ["complexity", "homology", "plasmid"]:
           
            if self.maindict[filt] == "unset":
                outline.append("-")
            elif readname in self.maindict[filt]:
                outline.append("1")
            else:
                outline.append("0")
        return "\t".join(outline)
                
                

def main():
    start_time_stamp = datetime.now()
    start_time = start_time_stamp.strftime("%Y/%m/%d, %H:%M:%S")
    
    sys.stdout.write('Started at: %s\n' % start_time)
    sys.stdout.flush()
    
    args = parser_io()
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # importing taxonomy
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    nodepath = os.path.join(args.taxonomy, "nodes.dmp")
    namepath = os.path.join(args.taxonomy, "names.dmp")
    merged_path = os.path.join(args.taxonomy, "merged.dmp")
    
    nodes = taxo.NodeTable(nodepath, merged_path)
    names = taxo.NameTable(namepath)
    
    names.namedict["0"] = taxo.Name("0", "unmapped", "unmapped", "unmapped")
    nodes.nodedict["0"] = taxo.Node("0", "1", *["unmapped"]*17)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # initiating filters
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    filters = Filter()
    
    if args.complexity:
        filters.add_filter(args.complexity, "complexity")
    if args.homology:
        filters.add_filter(args.homology, "homology")
    if args.plasmid:
        filters.add_filter(args.plasmid, "plasmid")
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # samfile parsing
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    samfile = open(args.insam)
    
    header = ["#flag", "name", "taxid", "length", "stats", "complexity", "homology", "plasmid"]
    
            
    readname = ""
    stat_line = ""
    cou = 0
    outfile = open(args.out, "w")
    outfile.write("\t".join(header) + "\n")
    
    for line in samfile:
        if line.startswith("@"):
            continue
        samline = SamLineXA(line)
        cou +=1
        
        #if cou >5:
        #    break
        
        if not readname:
            readname = samline.query_name
            stat_line = StatLine(readname)
            
        if samline.query_name != readname:
            # postprocessing stat_line
            stat_line.do_lca(nodes.nodedict)
            
            # applying filters
            filter_line = filters.filter_readname(stat_line.name)
            
            outfile.write(stat_line.__str__() + "\t" + filter_line + "\n")
            
            # creating new stat_line
            readname = samline.query_name
            
            stat_line = StatLine(readname)
            
        # fill into existing statline
        stat_line.add_samline(samline) 
        
    #catch the last element!
    stat_line.do_lca(nodes.nodedict)
    
    # applying filters
    filter_line = filters.filter_readname(stat_line.name)
    
    outfile.write(stat_line.__str__() + "\t"+ filter_line + "\n")
       
    samfile.close()

    outfile.close()    

    end_time_stamp = datetime.now()
    end_time = end_time_stamp.strftime("%Y/%m/%d, %H:%M:%S")
    sys.stdout.write('Finished at: %s\n' % end_time)
    sys.stdout.write('Duration: %s seconds\n'%(end_time_stamp - start_time_stamp))
    sys.stdout.flush()
    
    

if __name__ == "__main__":
    main()

        
        