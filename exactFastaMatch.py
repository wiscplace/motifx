#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
@Program: exactPeptideMatch.py

@Purpose: Given a list of peptides, gene names and a fasta file search for exact 
          matches and return gene name.
          
          The gene names will be used to limit the fasta search space, by keeping
          only those genes sequences in the fasta file that are in the gene list.
          
          This is an exact match search only!
          
          This has been tested only for peptide sequences but may also work
          for short DNA sequences.
          
@Input: a sequence file, gene name file, fasta file and output file name.

        sequence file is a text file (program written to take output of motifx.py)
        formatted as :
            
            PRARSSSVSNAAL,Induced,...R..S.S....
            LRERSRSNSSALA,Induced,...R..S.S....
            GSERTRSISFSKL,Induced,...R..S.S....
        
        but any text file with the peptide sequence in the first column should work.
        
        gene name file, this is the list of genes to search, this will be used
        to exclude genes from the fasta file. One gene name per line.
            YGL076C
            YNR047W
            YPL085W
        
        fasta file is a standard fasta file, but expected to be like the SGD
        orf file:  orf_trans_all_R64-2-1_20150113.fasta
        
        

@Dependencies: Python 3
               BioPython      
               
@Output: gene name with amino acid sequence and motif
@author: Mike Place
@Date:   12/02/2016
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict 
import argparse	
import os
import re
import sys

class exactFastaMatch(object):
    """ A sequence matching object.
    
    Attributes:
        genes:    Gene file to limit fasta search space.
        fastaRec: Fasta file to search for exact match.
        sequence: dictionary, key = 'STSSSPSSSPMSA', values = {'group': 'Induced', 'motif': '......S.SP...'}
        result:   matched gene name and sequence
    """

    def __init__(self, fasta, seq, gene):
        """ Contruct exactFastaMatch object """
        self.genes    = self.getGeneList(gene)      # load gene list 
        self.fastaRec = self.filterFasta(fasta)     # make a reduced list of seqIO objects (fasta sequences)
        self.sequence = self.getQuerySeq(seq)       #
        self.result   = []
    
    def getQuerySeq(self, seq):
        """ parse and return a dict of query sequences """
        seqlst = defaultdict(dict)
        
        with open(seq, 'r') as file:
            for ln in file:
                ln = ln.rstrip()
                row = ln.split(',')
                seqlst[row[0]]['group'] = row[1]
                seqlst[row[0]]['motif'] = row[2]
        return seqlst
    
    def getGeneList(self, gene):
        """ Read in gene list from text file """
        gnlst = []
        with open(gene,'r') as file:
            for g in file:
                g = g.rstrip()
                gnlst.append(g)
        return gnlst
    
    def filterFasta(self,fasta):
        """ select only those fasta sequences that match the names in the gene list """
        # open and create a list of all input fasta sequences
        record      = list(SeqIO.parse( fasta, 'fasta'))
        filteredRec = []
        # loop through all sequence records and match name to gene name
        # if names match keep record
        for rec in record:
            if rec.name in self.genes:
                filteredRec.append(rec)
        
        return filteredRec 
        
    def matchSeq(self):
        """ Attempt to match the sequences to the reduced gene list"""
        # loop through short sequence to match to fasta records
        for query in self.sequence.keys():
            for subject in self.fastaRec:            
                value = subject.seq.find(query)    # value is the index of the start position
                if value >= 0:
                    start = int(value)
                    end   = int(value) + 13
                    line = subject.name + '\t' + subject.seq[start:end] + '\t' + self.sequence[query]['motif'] + '\t' + self.sequence[query]['group']
                    self.result.append(line)
                    
    def writeResult(self, outfile):
        """ Write match gene names and motif information to file"""
        with open(outfile, 'w') as out:
            for r in self.result:
                out.write('%s\n' %(r))
        out.close
                     
def main():
    """
    Main, parse cmd line args and call writeRefFasta. 
    """
#******************************************************************************
# Command line args
#******************************************************************************    
    cmdparser = argparse.ArgumentParser(description="Match a short sequence to a gene in fasta file.",
                                        usage='%(prog)s -f <fasta file> -g <gene file> -s <sequence file>' ,prog='exactFastaMatch.py' )
    cmdparser.add_argument('-f', '--fasta',  action='store', dest="FASTA",  help='Required: fasta input file (.fasta, .fsa, .fa )', metavar='')
    cmdparser.add_argument('-d', '--detail', action='store_true', dest="DETAIL" , help='Print a more detailed description')
    cmdparser.add_argument('-g', '--gene',   action='store', dest="GENE", help='Required: Gene list file', metavar='')
    cmdparser.add_argument('-o', '--out',    action='store', dest='OUT',  help='Required: Output file name', metavar='')
    cmdparser.add_argument('-s', '--seq',    action='store', dest='SEQ',  help='Required: Query sequences file', metavar='')             
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)

#******************************************************************************
# Check and define input variables
#******************************************************************************    
    # Print detailed program information to screen
    if cmdResults['DETAIL']:
        print("\n\texactFastaMatch.py\n")
        print("\tPurpose: Exact match a short sequence to a gene in fasta file  ")
        print("\n\tRequired parameters:")
        print("\t  -f fasta file, expected to be an orf file.")
        print("\t\t where the gene name is right after the >")
        print("\t\t example:  >YAL001C ")
        print("\t  -g gene file, one gene name per line.")
        print("\t     Used to reduce the fasta search space.")
        print("\t\t one gene name per line")
        print("\t\t YAL001C")
        print("\t\t YAL002W")        
        print("\t  -o output file name")
        print("\t  -s sequence file, short sequences one per line.")        
        print("\t\t one short sequence per line, expect a CSV file")
        print("\t\t any extra columns will be ignored")
        print("\t\t PRARSSSVSNAAL")
        print("\t\t LRERSRSNSSALA")
        print("\t\t GSERTRSISFSKL")
        print("\n\toutput file:")
        print("\t Matched gene names in text file\n\n")
        sys.exit(1)
    
    # get fasta sequence file to query
    if cmdResults['FASTA'] is not None:
        fastaFile = cmdResults['FASTA']
        
    # get short sequences to find in the fasta sequences.
    if cmdResults['SEQ'] is not None:
        seqFile = cmdResults['SEQ']
        
    # get list of gene names to use to filter fasta sequence
    if cmdResults['GENE'] is not None:
        geneFile = cmdResults['GENE']   
    
    if cmdResults['OUT'] is not None:
        outfile = cmdResults['OUT']
    
    #for seq in SeqIO.parse(str(fastaFile), "fasta"):
    #    print(seq)
    data = exactFastaMatch(fastaFile, seqFile, geneFile)
    data.matchSeq()
    data.writeResult(outfile)
    
if __name__ == "__main__":
    main()

