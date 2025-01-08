#!/usr/bin/env python3

import os
import sys
import argparse
from datetime import datetime
import subprocess
import re
import numpy as np
import resource


import taxonomy_NCBI as taxo
import from_readstats_to_kmers as km

'''
Input:
* sample readstat file
* sample primary (no secondary/supplementary alignments, no unaligned) nohum (no human primary alignments) samfile sorted by reference

Output:
* taxrank
* taxid
* taxname
* count reads
* median readlen
* median alignment coverage
* median NM identity
* kmer index (num kmers/num kmers expected)

light variant of creae_taxstats_from_readstats.py
changed format of calculated statstics
deleted kraken stats, because other type of analysis. I could not calculate kmers and other future stats without mapping coordinates
'''

'''
ToDo:


do extra Table or two summarizing DB and filters
--mean,min,max of % assembly filter coverage in a taxon
--add number of kmers all and uniq per taxon (min,mean,max)
--add number genomes in this taxon (family, genus, species)



Done:
--align_len is now derived from readstat is align_len_onread, avoiding over 100% coverage
--NM is calcualted as: nm_percent = float(sum([x.nm for x in main_matches]))/float(align_len + delet)*100.

-- taxonomy and database fasta are now arguments, not hardcoded

--correct Kmer-Index: 
  -- do not add expected kmers as a single long alignment
  -- big issue: do so, that you can choose reads from a list of surviving when calculating the kmer index!!!

Neue Spalten:
--Reads_Survived_NM_and_Coverage
--Reads_Survived_Plasmid_Homology_Complexity
--Reads_Survived_All
--TaxonSurvived ?? –nur Reads, bei denen Kmer-Filter für Reads_Survived_All > Threshold)

'''


'''
 ~/my_scripts/Pipeline_v0.5/bin/from_readstats_to_taxstats_kmers.py SRR8289045_bwa_a_sorted_v05.readstat species SRR8289045_bwa_a_sorted_metagen_primary_byref.sam 32 SRR8289045_bwa_a_sorted_v05.taxstat
'''

def parser_io():
    parser = argparse.ArgumentParser(description="creating a taxstat file from readstats")
    parser.add_argument("readstats", help="readstats file")
    parser.add_argument("rank", help="taxonomic rank for the analysis")
    parser.add_argument("samfile", help="sorted (by position) sam file")
    parser.add_argument("k", type = int, help="k-mer size")
    parser.add_argument("outfile", help="output file")
    parser.add_argument("taxonomy", help="path to NCBI taxonomy file")
    parser.add_argument("database_fasta", help="path to the database fasta file")
    
    
    parser.add_argument("--nm_cutoff", help="nm cutoff")
    parser.add_argument("--coverage_cutoff", help="coverage cutoff")
    return parser.parse_args()



def log_time_ram(message):
    time_stamp = datetime.now()
    time_string = time_stamp.strftime("%Y/%m/%d, %H:%M:%S")
    
    ram = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) 
    ram_str = "%d kb"%ram
    if ram>1000 and ram<1000000:
        ram_str = "%d Mb"%(ram/1000.)
    elif ram>1000000:
        ram_str = "%d Gb"%(ram/1000000.)
    sys.stdout.write('%s: %s Peak RAM usage: %s\n'%(message,time_string, ram_str))
    sys.stdout.flush()
    return

#------------------------------------------------------------
def get_rank_taxid(taxid, rank, taxnodes_dict):
    node = taxnodes_dict[taxid]
    if node.taxid=="0":
        return "0"
    if node.rank == rank:
        return node.taxid
    if node.taxid=="1":
        return "-1" # no node at this rank
    return get_rank_taxid(node.parent, rank, taxnodes_dict)


def get_to_the_root(node_dict, node, way):
    way.append(node.taxid)
    if node.taxid=="1":
        way.reverse()
        return way
    return get_to_the_root(node_dict, node_dict[node.parent], way)

#------------------------------------------------------------

def make_stat(arr):
    if len(arr) == 7:
        return StatBWA(arr)
    else:
        sys.stderr.write("invalid statistics: %s\nExiting"%string)
        exit(0)

# both stat classes are used for single statistic and for statistic if a complete read or a read pair
class StatBWA:
    def __init__(self, arr):
        self.readlen, self.flag, self.taxid, self.nm, self.align_len, self.insert, self.delet = arr
        self.readlen = int(self.readlen)
        self.nm = int(self.nm)
        self.align_len = int(self.align_len)
        self.insert = int(self.insert)
        self.delet = int(self.delet)
        
    def as_array(self):
        return [self.readlen, self.nm, self.align_len, self.insert, self.delet]
    
    def __str__(self):
        strings = [str(x) for x in [self.flag, self.readlen ,self.nm, self.align_len, self.insert, self.delet]]
        return ":".join(strings)



class ClassLine:
    def __init__(self, line):
        self.arr = line.strip().split("\t")
        self.classified = self.arr[0]
        self.readname = self.arr[1]
        self.taxid = self.arr[2]
        self.length = self.arr[3].split("|")
        self.raw_stats = self.arr[4].split(" |:| ")
        self.complexity = int(self.arr[5])
        self.homology = int(self.arr[6])
        self.plasmid = int(self.arr[7])
        
        self.stats = []
        self.stat = '' # a single statistic of a read pair or a complete read 
        for fr_len, fr_stat in zip(self.length, self.raw_stats): # fr = forward reverse 
            arr = fr_stat.strip().split(" ")
            for st in arr:
                stat_arr = [fr_len] + st.strip().split(":")
                stat = make_stat(stat_arr)
                self.stats.append(stat)
                
    
    def calculate_stat_bwa(self):
        main_matches = [x for x in self.stats if x.flag=="M"]
        if  main_matches:
            readlen = sum([int(i) for i in self.length])
            #readlen = sum([x.readlen for x in main_matches])
            align_len = float(sum([x.align_len for x in main_matches]))
            align_len_percent = align_len/float(readlen)*100.
            insert = sum([x.insert for x in main_matches])
            delet = sum([x.delet for x in main_matches])
            nm_percent = float(sum([x.nm for x in main_matches]))/float(align_len + delet)*100. # add deletion to get alignment length on both read and reference
            self.stat = StatBWA([readlen, "M", self.taxid, nm_percent, align_len_percent, insert, delet]) # lca taxid!
        else:
            other_matches = [x for x in self.stats if x.flag!="U"]
            unmapped = [x for x in self.stats if x.flag=="U"]
            if other_matches:
                sys.stderr.write("secondary match without primary match!!!\n%s\nExiting"%"\t".join(self.arr))
                exit(0)
            if unmapped:
                self.stat = unmapped[0] 
            else:
                sys.stderr.write("no recognisable bwa matches!!!\n%s\nExiting"%"\t".join(self.arr))
                exit(0)     



# ToDo: 
# 1. create full ReadbasedStat dictionary
# 2. parse samfile with kmer beast
#      if the species in ReadbasedStat dictionary:
#            calculate the kmer index (rewrite the kmer beast in a simple "give me all sequences" thing, and then work with sequences in form key: readname, value: sequence covered in alignment)
class ReadbasedStats:
    def __init__(self, taxid):
        self.taxid = taxid
        self.readnames = [] # later to filter out surviving reads for kmer filter
        self.readlen = []
        self.align_coverage = []
        self.nm = []
        self.complexity = []
        self.homology = []
        self.plasmid = []
    
    def add_read(self, class_line_obj):
        clo = class_line_obj
        self.readnames.append(clo.readname)
        self.readlen.append(clo.stat.readlen)
        self.align_coverage.append(clo.stat.align_len)
        self.nm.append(clo.stat.nm)
        self.complexity.append(clo.complexity)
        self.homology.append(clo.homology)
        self.plasmid.append(clo.plasmid)
        
    def add_line(self, array):
        readname,readlen,align_cov,nm,compl,hom,plasm = array
        self.readnames.append(readname)
        self.readlen.append(readlen)
        self.align_coverage.append(align_cov)
        self.nm.append(nm)
        self.complexity.append(compl)
        self.homology.append(hom)
        self.plasmid.append(plasm)
        
        
    def subset_NC(self, nm_cutoff, cov_cutoff): # subsetting reads surviving nm and coverage filtering
        filtered = ReadbasedStats(self.taxid)
        for line in zip(self.readnames, self.readlen, self.align_coverage, self.nm, self.complexity, self.homology, self.plasmid):
            if line[3]<=nm_cutoff and line[2]>=cov_cutoff:
                filtered.add_line(line)
        return filtered
        
    def subset_HPC(self): # subsetting reads surviving homology, plasmid and complexity filters
        filtered = ReadbasedStats(self.taxid)
        for line in zip(self.readnames, self.readlen, self.align_coverage, self.nm, self.complexity, self.homology, self.plasmid):
            if line[4]==0 and line[5]==0 and line[6]==0:
                filtered.add_line(line)
        return filtered
    
    def as_array(self):
        count = len(self.readnames)
        if count>0:
            out = [count, np.median(self.readlen), np.median(self.align_coverage), np.median(self.nm), sum(self.complexity), sum(self.homology), sum(self.plasmid)]
        else:
            out = [count] + ["-"]*6
        return out
    
    def __str__(self):
        out = self.as_array()
        return "\t".join([str(x) for x in out])
        
#-----------------------------------------------------------------
        
def get_tax_route(ranks, taxid, nodedict, namedict):
    
    # get up to the next important rank
    node = nodedict[taxid]
    rank = node.rank
    while rank not in ranks:
        if node.taxid in ["0", "1"]:
            rank = "superkingdom"
            break
        node = nodedict[node.parent]
        rank = node.rank
        
    way_to_root = []
    # search names
    pos_rank = ranks.index(rank)
    
    for r in ranks[:pos_rank]:
        r_tid = get_rank_taxid(taxid, r, nodedict)
        name = "no_" + r
        if r_tid in namedict:
            name = namedict[r_tid].name
        way_to_root.append(name)    
    way_to_root = "|".join(way_to_root)
    return way_to_root


def output_line(taxid, nodedict, namedict, main_ranks):
    '''
    * taxrank
    * taxid
    * taxname
    * taxroute
    '''
    
    taxrank = nodedict[taxid].rank
    name = namedict[taxid].name
    route = get_tax_route(main_ranks, taxid, nodedict, namedict)
    out_arr = [taxrank, taxid, name, route]
    return out_arr



def main():
    start_time_stamp = datetime.now()
    start_time = start_time_stamp.strftime("%Y/%m/%d, %H:%M:%S")
    
    sys.stdout.write('Started at: %s\n' % start_time)
    sys.stdout.flush()
    
    args = parser_io()
    taxonomy=args.taxonomy
    db_fasta = args.database_fasta
    
    
    log_time_ram("Importing taxonomy")
    nodepath = os.path.join(taxonomy, "nodes.dmp")
    namepath = os.path.join(taxonomy, "names.dmp")
    merged_path = os.path.join(taxonomy, "merged.dmp")
    
    nodes = taxo.NodeTable(nodepath, merged_path)
    names = taxo.NameTable(namepath)
    
    
    # new to v0.4
    names.namedict["0"] = taxo.Name("0", "unclassified", "unclassified", "unclassified")
    nodes.nodedict["0"] = taxo.Node("0", "1", *["unclassified"]*17)
    
    names.namedict["-1"] = taxo.Name("-1", "no_rank", "no_rank", "no_rank")
    nodes.nodedict["-1"] = taxo.Node("-1", "1", *["no_rank"]*17)
    
    all_ranks = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"] # and MAX deep maybe later
    human_taxid = get_rank_taxid("9606", args.rank, nodes.nodedict)
    ignoring_taxids = [human_taxid, "0", "-1"]
    
    
    read_dict = dict() # key: readname, value: rank read taxid (after LCA given in readstat)
    #rank_dict = dict() # key: ranked taxid, value: list of stats
    rank_dict = dict() # key: ranked taxid, value: ReadbasedStats element
    cou = 0
    
    
    log_time_ram("Processing readstats")
    
    with open(args.readstats) as readstats:
        for line in readstats:
            if line.startswith("#"):
                continue
            
            #cou += 1
            #if cou > 100000:
            #    break
            
            l = ClassLine(line)
            rank_taxid = get_rank_taxid(l.taxid, args.rank, nodes.nodedict)
            
            l.calculate_stat_bwa()
                
            if not rank_taxid in rank_dict:
                rank_dict[rank_taxid] = ReadbasedStats(rank_taxid)
            
            rank_dict[rank_taxid].add_read(l)
    
    
    #target_readnames = [v.readnames for k,v in rank_dict.items() if k not in ignoring_taxids]
    #flat_target_readnames = [x for sublist in target_readnames for x in sublist]
    target_read_dict = dict() # very funky! key:readname, val:readname
    for key,val in rank_dict.items():
        if key not in ignoring_taxids:
            for read in val.readnames:
                target_read_dict[read] = read
    
    log_time_ram("Length all reads of interest: %d"%(len(target_read_dict)))
    
    log_time_ram("Producing output")
    
    header_arr = ["#rank", "taxid", "name", "taxroute"]
    header_part = ["count", "readlen_median", "coverage_median", "nm_median", "complexity", "homology", "plasmid", "kmer_index"]
    
    
    for i in ["all", "HPC", "NC", "filtered"]:
        header_arr += ["_".join([x,i]) for x in header_part]
    
    header = "\t".join(header_arr)
    
    outf = open(args.outfile, "w")
    outf.write(header + "\n")
    outf.flush()
    
    outlines = []
    
    cou = 0
    for key,val in rank_dict.items():
        cou += 1
        
        cutval_NC = val.subset_NC(10,70)
        cutval_HPC = val.subset_HPC()
        cutval_all = cutval_NC.subset_HPC()
        
        outline = output_line(key, nodes.nodedict, names.namedict, all_ranks)
        
        for dataset in [val,cutval_HPC,cutval_NC, cutval_all]:
            kmer_index = "-"
            
            log_time_ram("tax %s"%(dataset.taxid))
            outline.extend(dataset.as_array())
            outline.append(str(kmer_index))
        
        outlines.append(outline)
    
    outlines.sort(key=lambda x: int(x[4]), reverse=True)
    for l in outlines:
        outf.write("\t".join([str(x) for x in l]) + "\n")
    
    end_time_stamp = datetime.now()
    end_time = end_time_stamp.strftime("%Y/%m/%d, %H:%M:%S")
    sys.stdout.write('Finished at: %s\n' % end_time)
    sys.stdout.write('Duration: %s seconds\n'%(end_time_stamp-start_time_stamp))
    sys.stdout.flush()
    
    
if __name__ == "__main__":
    main()

