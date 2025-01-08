#!/usr/bin/env python3

import os
import sys
import argparse

def parser_io():
    parser = argparse.ArgumentParser(description="creating a taxstat file from readstats")
    parser.add_argument("taxstats_species", help="taxstat summary species level")
    parser.add_argument("taxstats_genus", help="taxstat summary genus level")
    parser.add_argument("taxstats_family", help="taxstat summary family level")
    parser.add_argument("--kmer_cutoff", type=float, default=0.7, help="kmer_index_cutoff")
    parser.add_argument("--count_cutoff", type=int, default=10, help="read count cutoff")
    parser.add_argument("--fraction_cutoff", type=float, default=1e-7, help="read fraction cutoff")
    parser.add_argument("out", help="output report")
    return parser.parse_args()

def forceint(instring):
    val = 0
    try:
        val = float(instring)
    except:
        pass
    return val

def get_all_reads(taxstat):
    with open(taxstat) as inf:
        abscount = 0
        header = inf.readline().strip("\n").split("\t")
        for line in inf:
            arr = line.strip("\n").split("\t")
            count = int(arr[4])
            abscount += count
        
    return abscount


class TaxstatLine:
    def __init__(self, linedict, all_reads):
        self.ld = linedict
        self.taxroute = self.ld["taxroute"]
        self.count = int(self.ld["count_filtered"])
        self.kidx = forceint(self.ld["kmer_index_filtered"])
        self.name = self.ld["name"]
        
        self.outarray = [self.ld[x] for x in ["taxid", "count_filtered", "readlen_median_filtered", "coverage_median_filtered", "nm_median_filtered", "kmer_index_filtered"]]
        fract = "{:.2E}".format(float(self.count)/float(all_reads))
        self.outarray.insert(1,fract) # myList.insert(0, 'zero')
        
        self.species = "-"
        self.genus = "-"
        self.family = "-"
        
        self.rank = self.ld["#rank"]
        rank = self.rank
        if rank == "species":
            self.family, self.genus = self.taxroute.rsplit("|",2)[-2:]
            self.species = self.name
        if rank == "genus":
            self.family = self.taxroute.rsplit("|",1)[-1]
            self.genus = self.name
        if rank == "family":
            self.family = self.name
        
    def __str__(self):
        if self.rank == "family":
            out = [self.name,"",""] + self.outarray
        if self.rank == "genus":
            out = ["",self.name,""] + self.outarray
        if self.rank == "species":
            out = ["","",self.name] + self.outarray
        return "\t".join(out)
        
def dummy_famline_dict(family_name):
    test_keys = ["#rank", "taxid", "name", "count_filtered", "readlen_median_filtered", "coverage_median_filtered", "nm_median_filtered", "kmer_index_filtered", "taxroute"]
    test_values = ["family", "-", family_name, "0", "-", "-", "-", "-", "-"]
    famline_dummy = {test_keys[i]: test_values[i] for i in range(len(test_keys))}
    return famline_dummy
             

def parse_taxstat(infile):
    lines = [] 
    all_reads = get_all_reads(infile)
    with open(infile) as inf:
        header = inf.readline().strip("\n").split("\t")
        for line in inf:
            arr = line.strip("\n").split("\t")
            taxroute = arr[3]
            taxid = arr[1]
            parent = taxroute.split("|")[-1]
            ld = dict()
            for i,j in zip(header,arr):
                ld[i] = j
            lines.append(TaxstatLine(ld, all_reads))
    return lines

def main():    
    args = parser_io()
    
    all_reads = get_all_reads(args.taxstats_species)
    
    vals_spec = parse_taxstat(args.taxstats_species)
    vals_genus = parse_taxstat(args.taxstats_genus)
    vals_fam = parse_taxstat(args.taxstats_family)
    
    # ToDO here! Give a clear cutoffs, even if they are stupid!!!
    #-------------------------------------------------------
    fraction_read_cutoff = float(all_reads)*args.fraction_cutoff
    count_cutoff = max(fraction_read_cutoff, args.count_cutoff)
    current_fraction = float(count_cutoff)/float(all_reads)
    # ToDo: fix the message, show the warning for cutoffs under 10 reads.
    message = "Nr reads: %d Read count cutoff: %d corresponds to read fraction %s"%(all_reads, count_cutoff, "{:e}".format(current_fraction))
    #-----------------------------------------------------------
    print(message)
    
    # filtering 
    filtered_fams = []
    for l in vals_fam:
        if l.count >= count_cutoff and l.kidx > args.kmer_cutoff:
           filtered_fams.append(l)
    
    filtered_genus = []
    for l in vals_genus:
        if l.count >= count_cutoff and l.kidx > args.kmer_cutoff:
           filtered_genus.append(l)
    
    filtered_specs = []
    for l in vals_spec:
        if l.count >= count_cutoff and l.kidx > args.kmer_cutoff:
           filtered_specs.append(l)
    
    
    # getting target families
    families_fromspec = set([x.family for x in filtered_specs])
    families_fromgen = set([x.family for x in filtered_genus])
    families_fromfam = set([x.family for x in filtered_fams])
    
    target_families = set().union(families_fromspec, families_fromgen, families_fromfam)
    
    #getting target genuses
    genuses_fromgen = set([x.genus for x in filtered_genus])    
    genuses_fromspec = set([x.genus for x in filtered_specs])
    
    target_genuses = set().union(genuses_fromgen, genuses_fromspec)
    
    #getting target species
    target_specs = set([x.species for x in filtered_specs])
    #target_specs = set([x.species for x in vals_spec if x.genus in target_genuses]) # reports even non-target specs for target genuses
    
    # grouping report lines
    groups = []
    
    human_family = [x for x in vals_fam if x.name=="Hominidae"]
    human_genus = [x for x in vals_genus if x.name=="Homo"]
    human_species = [x for x in vals_spec if x.name=="Homo sapiens"]
    
    if human_family and human_genus and human_species:
        group = [human_family[0], human_genus[0], human_species[0]]
        groups.append(group)
    
    for fam in target_families:
        group = []
        
        famlines = [x for x in vals_fam if x.family==fam]
        if not famlines:
            family_dummy = TaxstatLine(dummy_famline_dict(fam), 1)
            group.append(family_dummy)
        else: 
            group.append(famlines[0])
        
        genuslines = [x for x in vals_genus if x.family==fam and x.genus in target_genuses]
        
        for gline in genuslines:
            group.append(gline)
            speclines = [x for x in vals_spec if x.family==fam and x.genus==gline.genus and x.species in target_specs]
            
            for sline in speclines:
                group.append(sline)
        
        groups.append(group)
    
    groups.sort(key = lambda group: max([x.count for x in group]), reverse=True)
    
    # output
    header = "\t".join(["Family", "Genus", "Species", "Taxid", "Fraction", "Readcount", "Readlength", "Read Coverage", "NM", "Kmer Index"])
    with open(args.out, "w") as outf:
        outf.write(header + "\n")
        for group in groups:
            for line in group:
                outf.write(line.__str__() + "\n")
                
    
if __name__ == "__main__":
    main()