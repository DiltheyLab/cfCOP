#!/usr/bin/env python

# for taxonomy files from 17.12.2019

import sys
import os
import gzip
import argparse

def parser_io():
    parser = argparse.ArgumentParser(description="this script includes classes for parsing and changing NCBI taxonomy")

    parser.add_argument("nodefile", help="NCBI Taxonomy nodes.dmp file. Same version as namefile")
    
    parser_add_argument("namefile", help="NCBI Taxonomy names.dmp file. Same version as nodefile")
    
    return parser.parse_args()

'''
class Node:
    def __init__(self, taxid, parent, rank, embl_code, division_id, inherited_division, gencode_id, inherited_gc, mitochondrial_gencode, inherited_MGC, genbank_hidden, subtree_root, comments, *args):
        self.taxid = taxid #0
        self.parent = parent #1
        self.rank = rank #2
        self.embl_code = embl_code #3
        self.division_id = division_id #4
        self.inherited_division = inherited_division #5
        self.gencode_id = gencode_id #6
        self.inherited_gc = inherited_gc #7
        self.mitochondrial_gencode = mitochondrial_gencode #8
        self.inherited_MGC = inherited_MGC #9
        self.genbank_hidden = genbank_hidden #10
        self.subtree_root = subtree_root #11
        self.comments = comments #12
        self.rest_positions = args
        
    def __str__(self):
        #return "\t".join([" : ".join([i,j]) for i,j in self.__dict__.items()])
        return "\t".join(self.__dict__.values())
    
    def output_string(self):
        return "\t|\t".join(self.__dict__.values()) + "\t\n"
    
    def create_child(self, child_taxid):
        # pseudo strain taxid is composed as follows:
        # <pseudostrain_nr> 0000000 <original_taxid>
        copy_node = Node(*self.__dict__.values())
        copy_node.taxid = child_taxid
        copy_node.parent = self.taxid
        copy_node.rank = "pseudo_" + self.rank
        return copy_node
'''

class Node:
    def __init__(self, taxid, parent, rank, embl_code, division_id, inherited_division, gencode_id, inherited_gc, mitochondrial_gencode, inherited_MGC, genbank_hidden, subtree_root, *args):
        self.taxid = taxid #0
        self.parent = parent #1
        self.rank = rank #2
        self.embl_code = embl_code #3
        self.division_id = division_id #4
        self.inherited_division = inherited_division #5
        self.gencode_id = gencode_id #6
        self.inherited_gc = inherited_gc #7
        self.mitochondrial_gencode = mitochondrial_gencode #8
        self.inherited_MGC = inherited_MGC #9
        self.genbank_hidden = genbank_hidden #10
        self.subtree_root = subtree_root #11
        self.rest_positions = args
        
    def __str__(self):
        #return "\t".join([" : ".join([i,j]) for i,j in self.__dict__.items()])
        return "\t".join(self.__dict__.values())
    
    def output_string(self):
        return "\t|\t".join(self.__dict__.values()) + "\t\n"
    
    def create_child(self, child_taxid):
        # pseudo strain taxid is composed as follows:
        # <pseudostrain_nr> 0000000 <original_taxid>
        copy_node = Node(*self.__dict__.values())
        copy_node.taxid = child_taxid
        copy_node.parent = self.taxid
        copy_node.rank = "pseudo_" + self.rank
        return copy_node

class MergedNodes:
    def __init__(self, filepath):
        self.nodelink = dict()
        with open(filepath) as infile:
            for line in infile:
                line = line.rstrip("\t|\n")
                from_node, to_node = line.split("\t|\t")
                self.nodelink[from_node] = to_node
        
        
class NodeTable:
    def __init__(self, path, merged):
        self.path = path
        self.merged = merged
        self.nodelist = []
        self.nodedict = dict()
        with open(self.path) as infile:
            for line in infile:
                line = line.rstrip("\t|\n")
                linearr = line.split("\t|\t")
                try:
                    node = Node(*linearr)
                except:
                    try:
                        new_id = merged.nodelink[linearr[0]]
                        node_args = [new_id] + linearr[1:]
                        node = Node(*node_args)
                        
                    except:
                        sys.stderr.write("Error: invalid taxonomic node: %s\nCheck the input taxonomic file nodes.dmp.\nLength of an input array: %d\nExiting...\n"%(" ".join(linearr), len(linearr)))
                        exit(0)
                    
                self.nodelist.append(node)
                self.nodedict[node.taxid] = node
                
    def add_node(self, node):
        if not node.taxid in self.nodedict:
            self.nodelist.append(node) 
            self.nodedict[node.taxid] = node

class Name:
    def __init__(self, taxid, name, unique_name, name_class):
        self.taxid = taxid
        self.name = name
        self.unique_name = unique_name
        self.name_class = name_class
    
    def __str__(self):
        return "\t".join(self.__dict__.values())
    
    def output_string(self):
        return "\t|\t".join(self.__dict__.values()) + "\t\n"
    
    def create_child(self, taxid):
        copy_name = Name(*self.__dict__.values())
        copy_name.taxid = taxid
        copy_name.name = "Pseudotaxonomic unit " + self.name
        copy_name.unique_name = ""
        return copy_name
        
        

class NameTable:
    def __init__(self, path):
        self.path = path
        self.namelist = []
        self.namedict = dict()
        with open(self.path) as infile:
            for line in infile:
                line = line.rstrip("\t|\n")
                linearr = line.split("\t|\t")
                try:
                    name = Name(*linearr)
                    self.namelist.append(name)
                    if not name.taxid in self.namedict:
                        self.namedict[name.taxid] = name
                    elif name.name_class == "scientific name":
                        self.namedict[name.taxid] = name
                except:
                    print(len(linearr), "should be 4")
                    sys.stderr.write("Error: invalid taxonomic name: " + " ".join(linearr) + "\nCheck the input taxonomic file names.dmp.\nExiting...\n")
                    exit(0)
    
    def add_name(self, name):
        if not name.taxid in self.namedict:
            self.namelist.append(name)
            self.namedict[name.taxid] = name


def open_compressed(infilepath):
    if infilepath.endswith('.gz'):
        f = gzip.open(infilepath, 'rt') # mode 'rt' opens file in an text format, default is binary!!!
        return f
    else:
        f = open(infilepath)
        return f


class Accession2TaxidLine:
    def __init__(self, accession, acc_version, taxid, gi):
        self.accession = accession
        self.acc_version = acc_version
        self.taxid = taxid 
        self.gi = gi
    
    def __str__(self):
        return "\t".join(self.__dict__.values())

def generate_accession2taxids(path):
    with open_compressed(path) as infile:
        infile.readline()
        for line in infile:
            linearr = line.strip().split("\t")
            try:
                acc2taxid = Accession2TaxidLine(*linearr)
                yield acc2taxid
            except:
                sys.stderr.write("Error: invalid line in %s. Please provide valid file\nLine: %s\nExiting"%(path," ".join(linearr)))
                exit(0)
            
        
def main():
    args = parser_io()
    a = NodeTable(args.nodefile)
    print(len(a.nodelist))
    for i in a.nodelist[:10]:
        print(i)
    all_avail_ranks = list(set([i.rank for i in a.nodelist]))
    print(all_avail_ranks)
    
    
 
if __name__ == "__main__":
    main()