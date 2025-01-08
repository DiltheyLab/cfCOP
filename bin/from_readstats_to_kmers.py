#!/usr/bin/env python3

import os
import sys
import argparse
from datetime import datetime
import subprocess
import re
import math

'''
This script calculates stats for one sample!
1. takes a readstat file (read statistics after LCA)
2. brings each read/read pair to a given rank
3. flags interesting reads (which have this rank)
4. parses the small sam file with only non-human main matches (flags 0 or 16) sorted by reference!!!!
5. for each read: is it an interesting read? does it have a flag?
6. if interesting read, count its kmers and theoretical on ref length to the certain taxid


As mainscript:
makes a table 
taxid, taxname, read_count, num kmers, expected kmers, kmer index

as module:
work on this!!! to integrate into taxstats

'''
sys.path.append('/home/tyshaiev/my_scripts/refseq_db_download_and_indexing')
import taxonomy_NCBI as taxo

def parser_io():
    parser = argparse.ArgumentParser(description="counting k-mers in a mapped sequence")
    parser.add_argument("readstats", help="readstats file")
    parser.add_argument("rank", help="taxonomic rank for the analysis")
    parser.add_argument("fasta", help="big database fasta")
    parser.add_argument("samfile", help="sorted (by position) sam file")
    parser.add_argument("k", type = int, help="k-mer size")
    parser.add_argument("outfile", help="output file")
    
    return parser.parse_args()

class ReadstatLine:
    def __init__(self, line):
        linearr = line.strip("\n").split("\t")
        self.readname = linearr[1]
        self.lca_taxid = linearr[2]
        
class SamLine:
    def __init__(self, line):
        linearr = line.strip("\n").split("\t")
        self.name = linearr[0]
        self.flag = linearr[1]
        self.refname = linearr[2]
        self.pos = int(linearr[3])
        self.cigar = CIGAR(linearr[5])
        self.on_ref = self.cigar.len_on_ref
        #self.read_seq = linearr[9]
        self.start_on_ref = self.pos - 1
        '''
        self.start_on_read = 0
        for pos, match in self.cigar.tuples: # check this!!!
            if match == 'S':
                self.start_on_ref += pos
                self.start_on_read += pos
            elif match in 'MDN=X':
                break
        '''
        self.stop_on_ref = self.start_on_ref + self.on_ref - 1
        #self.cutread = self.read_seq[self.start_on_read:self.start_on_read + self.cigar.len_on_read]
        
class CIGAR:
    def __init__(self, string):
        pattern = re.compile('([MIDNSHPX=])')
        vals = pattern.split(string)[:-1]
        self.tuples = list(([int(vals[i]), vals[i + 1]] for i in range(0, len(vals), 2)))
        self.len_on_ref = self.summarize('MDN=X') # all that consumes reference
        self.len_on_read = self.summarize('MI=X') # all that consumes query but S
        
    def summarize(self, patterns):
        arr = [int(pos) for pos, match in self.tuples if match in list(patterns)]
        return sum(arr) # example: list('MDN=X')
'''
class Sequence:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
'''            
def reverse_complement(seq):
    revdict = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N", "W":"W", "S":"S", "M":"K", "K":"M", "R":"Y", "Y":"R", "B":"V", "V":"B", "D":"H", "H":"D"}
    rev_compl = ""
    reversed_seq = seq[::-1]
    for char in reversed_seq:
        rev_compl += revdict[char]
    return rev_compl
    
def log_time(message):
    time_stamp = datetime.now()
    time_string = time_stamp.strftime("%Y/%m/%d, %H:%M:%S")
    sys.stdout.write('%s: %s\n'%(message,time_string))
    sys.stdout.flush()
    return

def get_rank_taxid(taxid, rank, taxnodes_dict):
    node = taxnodes_dict[taxid]
    if node.taxid=="0":
        return "0"
    if node.rank == rank:
        return node.taxid
    if node.taxid=="1":
        return "-1" # no node at this rank
    return get_rank_taxid(node.parent, rank, taxnodes_dict)

def get_ref_sequence(fasta, refname): # fasta should be indexed by SamTools!
    bashline = ['samtools', 'faidx', fasta, refname]
    process = subprocess.Popen(bashline, stdout=subprocess.PIPE)
    refseq = str(process.stdout.read(), 'utf-8')
    sequence_list = refseq.strip().split("\n")[1:]
    sequence = ""
    for i in sequence_list:
        sequence += i
    return sequence

# will be deprecated
def get_kmers(sequence, k):
    # returns set of kmers
    kmers = set()
    for i in range(0, len(sequence)-k+1):
        kmer = sequence[i:i+k]
        kmer = kmer.upper()
        rc_kmer = reverse_complement(kmer)
        if kmer < rc_kmer:
            kmers.add(kmer)
        else:
            kmers.add(rc_kmer)
    return kmers

class ReadBeast:
    def __init__(self): # lookup dict is the dict with key: readname and val: ranked taxid. Produced from readstats
        # body (stuffed by digestion)
        self.rank_reads_uniq = dict() # key: taxid, val: number reads
        self.rank_reads_all = dict()
        # general
        self.empty = 1 # is th KmerBeast empty
        # belly (emptied or changed by digestion)
        self.ref_reads_uniq = set()
        self.ref_reads_all = 0
        self.taxid = ""
    
    def eat_samline(self, samline, newref = "", taxid = ""): #if you feed both seq and ref, it will digest and initiate 
        sl = samline
        
        if newref: # feeder will always feed the reference suitable for the read
            if not self.empty:
                self.digest()
            else:
                self.empty = 0
            self.initiate(sl, newref, taxid)
        self.ref_reads_uniq.add(sl.name)
        self.ref_reads_all += 1
    
    def digest(self):
        if not self.taxid in self.rank_reads_uniq:
            self.rank_reads_uniq[self.taxid] = set()
        self.rank_reads_uniq[self.taxid] = self.rank_reads_uniq[self.taxid].union(self.ref_reads_uniq)
        if not self.taxid in self.rank_reads_all:
            self.rank_reads_all[self.taxid] = 0
        self.rank_reads_all[self.taxid] += self.ref_reads_all
        
    
    def initiate(self, sl, newref, taxid): # veeeery important is that samfile sorted byref
        self.ref_reads_uniq = set()
        self.ref_reads_all = 0
        self.taxid = taxid
        
    

class KmerBeast:
    def __init__(self, k): # lookup dict is the dict with key: readname and val: ranked taxid. Produced from readstats
        # body (stuffed by digestion)
        self.rank_kmer_dict = dict() 
        # general
        self.empty = 1 # is th KmerBeast empty
        self.k = k
        # belly (emptied or changed by digestion)
        self.ref = ""
        self.readname = ""
        self.taxid = ""
        self.ref_kmer_dict = dict() 
         
        
    def eat_samline(self, samline, newref = "", taxid = ""): #if you feed both seq and ref, it will digest and initiate 
        sl = samline
        
        if newref: # feeder will always feed the reference suitable for the read
            if not self.empty:
                self.digest()
            else:
                self.empty = 0
            self.initiate(sl, newref, taxid)
        
        aligned_ref = self.ref.seq[sl.start_on_ref:sl.stop_on_ref]
        self.fill_kmers(aligned_ref)
        
    # NOT FORGET maximal expected kmer num is sum over kmer dict values    
    def fill_kmers(self, sequence):
        k = self.k
        for i in range(0, len(sequence)-k+1):
            kmer = sequence[i:i+k]
            kmer = kmer.upper()
            rc_kmer = reverse_complement(kmer)
            kmer_id = min(kmer, rc_kmer)
            if not kmer in self.ref_kmer_dict:
                self.ref_kmer_dict[kmer_id] = 0
            self.ref_kmer_dict[kmer_id] += 1
        
        
    def initiate(self, sl, newref, taxid): # veeeery important is that samfile sorted byref
        self.ref = newref
        self.taxid = taxid
        self.ref_kmer_dict = dict()
        

    # NOT FORGET to make last digest after forloop in main
    def digest(self):
        # fill rank_kmer_dict
        #log_time("Digesting %s..."%self.taxid)
        if not self.taxid in self.rank_kmer_dict:
            self.rank_kmer_dict[self.taxid] = dict()
        for kmer, count in self.ref_kmer_dict.items():
            if not kmer in self.rank_kmer_dict[self.taxid]:
                self.rank_kmer_dict[self.taxid][kmer] = 0
            self.rank_kmer_dict[self.taxid][kmer] += count
        # --------- DEBUG ---------
                
                

                
class RipleyBeast:
    # ToDO: 
    # get first, last start, difference between them, contig len
    def __init__(self):
        # body (stuffed by digestion)
        # simplified RipleyK value is summ/reflen
        self.reflen_dict = dict()
        self.summ_dict = dict()
        # general
        self.empty = 1
        #self.search_radius = search_radius 
        # belly (emptied or changed by digestion)
        self.starts = []
        self.reflen = 0
        self.taxid = ""
    
    def eat_samline(self, samline, newref = "", taxid = ""):
        sl = samline
        
        if newref:
            if not self.empty:
                self.digest()
            else:
                self.empty = 0
            self.initiate(sl, newref, taxid)
        self.starts.append(samline.pos)
        
    def digest(self):
        # here digest is calculate the one-dimensional ripley K index
        summ = 0
        if len(self.starts)>1:
            
            search_radius = self.reflen/(len(self.starts)-1)
            
            for pos in range(len(self.starts)-1):
                x = self.starts[pos]
                y = self.starts[pos+1]
                neighbor_dist = abs(x-y)
                if neighbor_dist < search_radius:
                    summ += abs(x-y)
        else:
            self.reflen = 0 # here: dont count reference if only one read!
        
        if not self.taxid in self.reflen_dict:
            self.reflen_dict[self.taxid] = 0
        self.reflen_dict[self.taxid] += self.reflen
        if not self.taxid in self.summ_dict:
            self.summ_dict[self.taxid] = 0
        self.summ_dict[self.taxid] += summ
        
    def initiate(self, samline, newref, taxid):
        self.starts = []
        self.reflen = len(newref.seq)
        self.taxid = taxid
        
    
class CoverBeast:
    def __init__(self):
        pass
'''
1. Shannon index needs:
    * kmer dict
    
2. Ripley K needs:
    * starts of all reads

3. coverage stat needs:
    * list of positions which are taken
    * summary length of all reads
    
    Idea: make clases with body and belly. Belly takes a samfile line, if it is processed completely, metric goes to the body. Body is a dict where keys are taxids. 
    
    
'''

def test_soft_clips(samline, ref):
    print("-----------------------")
    #print("Read:", samline.read_seq)
    #print("Ref:", ref.seq[samline.pos:samline.pos+len(samline.read_seq)])
    print("Ref_cut:\t", ref.seq[samline.pos-1:samline.pos + samline.on_ref-1])
    print("Read_cut:\t", samline.cutread)
    print("CIGAR", samline.cigar.tuples, samline.pos, samline.flag)
    return


'''
# generator (SIMPLE), yields reference, readname, sequence covered on reference
def parse_sam_simple(insam, fasta):
    current_ref = Sequence("","")
    with open(insam) as samfile:
        for line in samfile:
            if line.startswith("@"):
                continue
                
            sl = SamLine(line)
            if sl.refname == "*":
                continue
                
            if sl.refname != current_ref.name:
                # get new reference sequence
                seq = get_ref_sequence(fasta, sl.refname)
                #set new current ref
                current_ref = Sequence(sl.refname, seq)
            
            pos_left = sl.pos-1 # 1-based leftmost position translated to pythonic 0-based
            pos_right = pos_left + sl.on_ref
            covered_ref_part = current_ref.seq[pos_left:pos_right]
                
            yield sl.refname, sl.name, covered_ref_part
'''

# generator, yields reference, read count, number_kmers, expected_number_kmers for a certain contig
# interesting reads is a list with interesting readname, only reads from this list will be reported
# CHANGED in v0.4: kmers is not a set, but a dict 
def parse_samfile(insam, fasta, interesting_reads, k):
    current_ref = Sequence("","")
    ref_kmers = dict()
    ref_on_ref_len = 0
    read_count = 0
    first_read = ""
    with open(insam) as samfile:
        for line in samfile:
            if line.startswith("@"):
                continue
            
            sl = SamLine(line)
            if not sl.name in interesting_reads:
                continue
            if sl.refname == "*":
                continue
            
            
            if sl.refname != current_ref.name:
                
                if current_ref.name:
                    #harvest old ref
                    expected_kmers = max(0, ref_on_ref_len - k + 1)
                    yield(current_ref.name, first_read, ref_kmers, expected_kmers, read_count) # here you report first read name to be able to find read rank taxid in the dictionary
                    
                    #ref_kmers = set()
                    ref_kmers = dict()
                    ref_on_ref_len = 0
                    read_count =  0
                    
                #set first read of the reference
                first_read = sl.name
                #find new seq
                seq = get_ref_sequence(fasta, sl.refname)
                #set new current ref
                current_ref = Sequence(sl.refname, seq)
                
            #------------------
            #softclip test
            #if sl.start_on_ref > sl.pos:
                #print(sl.name)
            #test_soft_clips(sl, current_ref)
            #------------------
                
            pos_left = sl.pos-1 # 1-based leftmost position translated to pythonic 0-based
            pos_right = pos_left + sl.on_ref
            covered_ref_part = current_ref.seq[pos_left:pos_right]
            kmers = get_kmers(covered_ref_part, k)
            for kmer in kmers:
                if not kmer in ref_kmers:
                    ref_kmers[kmer] = 0
                ref_kmers[kmer] += 1
            ref_on_ref_len += sl.on_ref
            read_count += 1
            
        #catch last element
        expected_kmers = max(0, ref_on_ref_len - k + 1)
        if current_ref.name:
            yield(current_ref.name, first_read, ref_kmers, expected_kmers, read_count) # here you report first read name to be able to find read rank taxid in the dictionary

def calc_shannon_index(kmer_dict): # this kmer dict should be key: kmer, value: number of kmer occurances
    num_kmers = sum(kmer_dict.values())
    if num_kmers < 2:
        return 0
    #print(canonical_kmers)
    ps = [float(i)/float(num_kmers) for i in kmer_dict.values()]
    #print(ps)
    p_log_ps = [i*math.log(i) for i in ps]
    #print(p_log_ps)
    H = -sum(p_log_ps)
    return H/math.log(num_kmers)
    

#------------HERE START!!!!!!!!--------------------#
def open_compressed(infilepath):
    if infilepath.endswith('.gz'):
        f = gzip.open(infilepath, 'rt') # mode 'rt' opens file in text format, default is binary!!!
        return f
    else:
        f = open(infilepath)
        return f

class Sequence:
    def __init__(self, name, seq):
        self.name = name
        self.shortname = name.split(" ")[0]
        try:
            self.taxid = self.shortname.strip("sp").split("as")[0]
        except:
            self.taxid = "0"
        self.seq = seq
        self.length = len(self.seq)
        self.line_split = 80 # line split of NCBI
        
    def __str__(self):
        outlines = [">" + self.name]
        outlines.extend([self.seq[i:i+self.line_split] for i in range(0, len(self.seq), self.line_split)])
        return "\n".join(outlines) + "\n"


def fasta_read(infile):
    f = open_compressed(infile)
    if not f:
        sys.stdout.write('Infile should be fna or fna.gz. Please check your input!\nExiting\n')
        exit(1)
    STATE_NAME = 0
    STATE_SEQ = 1
    state = STATE_NAME
    name = ""
    seq = ""
    for line in f:
        if state == STATE_SEQ:
            if line.strip().startswith(">"):
                state = STATE_NAME
                yield(Sequence(name, seq))
            else:
                seq += line.strip()
        if state == STATE_NAME:
            name = line.strip()[1:]
            state = STATE_SEQ
            seq = ""
    yield Sequence(name, seq) # to catch the last element
    f.close()


# generator (SIMPLE), yields reference, readname, sequence covered on reference
def parse_sam_simple(insam, fasta_dict):
    current_ref = Sequence("","")
    with open(insam) as samfile:
        for line in samfile:
            if line.startswith("@"):
                continue
                
            sl = SamLine(line)
            if sl.refname == "*":
                continue
                
            if sl.refname != current_ref.name:
                if not sl.refname in fasta_dict:
                    continue
                # get new reference sequence
                seq = fasta_dict[sl.refname]
                #set new current ref
                current_ref = Sequence(sl.refname, seq)
            
            pos_left = sl.pos-1 # 1-based leftmost position translated to pythonic 0-based
            pos_right = pos_left + sl.on_ref
            covered_ref_part = current_ref.seq[pos_left:pos_right]
                
            yield sl.refname, sl.name, covered_ref_part


def parse_sam_simple_target(insam, fasta_dict, target_reads): # target reads is a dictionary!!! key:readname, val:readname
    current_ref = Sequence("","")
    with open(insam) as samfile:
        for line in samfile:
            if line.startswith("@"):
                continue
            
            readname = line.split("\t",1)[0]
            readname = target_reads.pop(readname, None)
            
            if not readname:
                continue
            
            sl = SamLine(line)
            
            '''
            if sl.refname == "*":
                continue
            '''
            
            if sl.refname != current_ref.name:
                if not sl.refname in fasta_dict:
                    continue
                # get new reference sequence
                seq = fasta_dict[sl.refname]
                #set new current ref
                current_ref = Sequence(sl.refname, seq)
            
            pos_left = sl.pos-1 # 1-based leftmost position translated to pythonic 0-based
            pos_right = pos_left + sl.on_ref
            covered_ref_part = current_ref.seq[pos_left:pos_right]
                
            yield sl.refname, sl.name, covered_ref_part


#------------HERE END!!!!!!!!--------------------#

    
    
    
    
def main():
    # constants-------------------------
    ignoring_taxids = ["9606", "0", "-1"]
    taxonomy="/gpfs/project/dilthey/projects/HematoDETECT/refseq_db_2019_12_26/taxonomy"
    
    #-----------------------------------
    
    args = parser_io()
    k = args.k
    
    start_time_stamp = datetime.now()
    log_time("Started at")
    
    # importing taxonomy
    log_time("Importing taxonomy")
    nodepath = os.path.join(taxonomy, "nodes.dmp")
    namepath = os.path.join(taxonomy, "names.dmp")
    merged_path = os.path.join(taxonomy, "merged.dmp")
    
    nodes = taxo.NodeTable(nodepath, merged_path)
    names = taxo.NameTable(namepath)
    
    
    names.namedict["0"] = taxo.Name("0", "unmapped", "unmapped", "unmapped")
    nodes.nodedict["0"] = taxo.Node("0", "1", *["unmapped"]*17)
    
    # magic
    # flag interesting reads from readstats
    log_time("Parsing readstats")
    interesting_reads = dict()
    with open(args.readstats) as statfile:
        for line in statfile:
            if line.startswith("#"):
                continue
            rstat = ReadstatLine(line)
            rank_taxid = get_rank_taxid(rstat.lca_taxid, args.rank, nodes.nodedict)
            if rank_taxid in ignoring_taxids:
                continue
            interesting_reads[rstat.readname] = rank_taxid
    
    
    log_time("Readstats: %d interesting reads"%len(interesting_reads))
    
    log_time("Beast test")
    kmer_beast = KmerBeast(args.k)
    read_beast = ReadBeast()
    ripley_beast = RipleyBeast()
    
    ref = Sequence("","")
    with open(args.samfile) as samfile:
        for line in samfile:
            if line.startswith("@"):
                continue
            
            sl = SamLine(line)
            if not sl.name in interesting_reads:
                continue
                
            if sl.refname == "*":
                continue
            
            if sl.refname != ref.name:
                refseq = get_ref_sequence(args.fasta, sl.refname)
                ref = Sequence(sl.refname, refseq)
                rank_taxid = interesting_reads[sl.name]
                kmer_beast.eat_samline(sl, ref, rank_taxid)
                read_beast.eat_samline(sl, ref, rank_taxid)
                ripley_beast.eat_samline(sl, ref, rank_taxid)
            
            else:
                kmer_beast.eat_samline(sl)
                read_beast.eat_samline(sl)
                ripley_beast.eat_samline(sl)
    
    kmer_beast.digest()
    read_beast.digest()
    ripley_beast.digest()
    
    #log_time("All reads eaten by KmerBeast: %d"%(kmer_beast.reads))
    #log_time("Uniq reads eaten by KmerBeast: %d"%(len(kmer_beast.uniq_reads)))
    
    log_time("Beast output")
    outfile = open(args.outfile, "w")
    
    outarr = ["#rank","taxid", "name", "num_reads_uniq", "num_reads_all", "uniq_kmers", "all_kmers", "kind", "sind", "Krip"]
    outfile.write("\t".join(outarr) + "\n")
    
    for taxid, kmers in kmer_beast.rank_kmer_dict.items():
        rank = nodes.nodedict[taxid].rank
        taxname = names.namedict[taxid].name
        num_reads_uniq = len(read_beast.rank_reads_uniq[taxid])
        num_reads_all = read_beast.rank_reads_all[taxid]
        uniq_kmers = len(kmers.keys())
        all_kmers = sum(kmers.values())
        kind = 0
        sind = calc_shannon_index(kmers)
        if all_kmers:
            kind = float(uniq_kmers)/float(all_kmers)
        Krip = "None"
        ripley_reflen =  ripley_beast.reflen_dict[taxid]
        if ripley_reflen:
            ripley_summ = ripley_beast.summ_dict[taxid]
            Krip = float(ripley_summ)/float(ripley_reflen)
        outarr = [rank, taxid, taxname, num_reads_uniq, num_reads_all, uniq_kmers, all_kmers, kind, sind, Krip]
        outfile.write("\t".join([str(x) for x in outarr]) + "\n")
        
    
    
    
    exit(0)
    
    log_time("Started parsing samfile")
    # parse small samfile
    kmer_dict = {x:dict() for x in interesting_reads.values()}
    expected_kmers = {x:0 for x in interesting_reads.values()}
    counts = {x:0 for x in interesting_reads.values()}
    
    for ref, first_read, kmers, exp_kmers, read_count in parse_samfile(args.samfile, args.fasta, interesting_reads.keys(), 32):
        rank_taxid = interesting_reads[first_read]
        #print(ref,rank_taxid, len(kmers), expected_kmers, read_count)
        #kmer_dict[rank_taxid] = kmer_dict[rank_taxid].union(kmers)
        # combining dicts here
        for kmer,num in kmers.items():
            if not kmer in kmer_dict[rank_taxid]:
                kmer_dict[rank_taxid][kmer] = 0
            kmer_dict[rank_taxid][kmer] += num
        expected_kmers[rank_taxid] += exp_kmers
        counts[rank_taxid] += read_count
        
        
    
    # output
    out_matrix = []
    header = ["#taxid","taxname","count","kmers","expected_kmers", "kmer_index", "shannon_index"]
    for taxid in kmer_dict.keys():
        taxname = names.namedict[taxid].name
        count = counts[taxid]
        kmers = len(kmer_dict[taxid])
        exp_kmers = expected_kmers[taxid]
        if exp_kmers>0:
            kmer_index = "{:.2f}".format(float(kmers)/float(exp_kmers))
        else:
            kmer_index = "0.00"
        shannon_index = "{:.2f}".format(calc_shannon_index(kmer_dict[taxid]))
        outline = [taxid, taxname, count, kmers, exp_kmers, kmer_index, shannon_index]
        out_matrix.append(outline)
    
    # writing output
    with open(args.outfile, "w") as outf:
        outf.write("\t".join(header) + "\n")
        for line in out_matrix:
            outf.write("\t".join([str(x) for x in line]) + "\n")
    
    
    end_time_stamp = datetime.now()
    end_time = end_time_stamp.strftime("%Y/%m/%d, %H:%M:%S")
    sys.stdout.write('Finished at: %s\n' % end_time)
    sys.stdout.write('Duration: %s seconds\n'%(end_time_stamp-start_time_stamp))
    sys.stdout.flush()
    
    
if __name__ == "__main__":
    main()
