#!/usr/bin/python

__author__ = "Andriy Sheremet"

from Bio import SeqIO, SeqRecord, Seq
from Bio.Alphabet import IUPAC
import getopt
import sys
import operator
import itertools as it
try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO


pars ={}

pars["infile"] = ""


def usage():
    print "\nThis is the usage function:\n"

    

def main(argv):
    
    #default parameters
    global pars
             
    try:                                
        opts, args = getopt.getopt(argv, "i:h", ["input=", "help"])
    except getopt.GetoptError:          
        usage()                         
        sys.exit(2)                     
    for opt, arg in opts:                
        if opt in ("-h", "--help"):      
            usage()
            supported_operations()
            sys.exit()                  
        elif opt in ("-i", "--input"):
            if arg:
                pars["infile"] = arg
            #print "Input file", arg                  
    #print("hehe")
    #print(pars['infile'])
    extract_features(pars)

    
def extract_features(params):
    #print("extract_features() is called")
    

    
    if params["infile"]:
        try:
            #
            handle = open(params["infile"], "rU")
        except:
            print "\nERROR: Input file doesn't exist"
            usage()
            sys.exit(1)
    else:
        try:
            #
            handle = sys.stdin
        except:
            print "\nERROR: Input file doesn't exist"
            usage()
            sys.exit(1)
            

    
    OUT_FILE = params["infile"]+".faa"
    #print("hehe")
    #print(type(faa_filename))
    #output_handle = open(faa_filename, "w")
    #print("haha")
    entries=list()
    headers=list()
    
    
    for key,group in it.groupby(handle,lambda line: line.startswith('LOCUS')):
            if not key:
                entries.append(list(group))
            else:
                headers.append(list(group))
                
    combined=map(lambda x1, x2: x1+x2, headers, entries)
    

    
    
    
    seqs_for_fasta = []
    
    for item in combined:
        try:
            f=StringIO("".join(item))
            record=SeqIO.read(f, 'genbank')
        except:
            print("Error: Could not parse genbank file ["+params["infile"]+"]")
            
            
        try:
            org = record.annotations['organism']
        except:
            org = ''

            #acc = record.annotations['accessions'][0]  ## -not needed for now
            #tax_line = ("; ").join(record.annotations['taxonomy'] ## -not needed for now, need to be tested if works
    
        for feature in record.features:
            # Each Genbank record may have several features, the program
            # will walk over all of them.
            qualifier = feature.qualifiers
        #     print("<<")
        #     print(qualifier)
        #     print(">>")
            if 'product' in qualifier and 'translation' in qualifier:#we've got a protein
                
                id_ = qualifier['protein_id'][0]
                desc = qualifier['product'][0]
                sq = Seq.Seq(qualifier['translation'][0], IUPAC.protein)
                if org:
                    srec = SeqRecord.SeqRecord(sq, id=id_, description=desc+" ["+org+"]")
                else:
                    srec = SeqRecord.SeqRecord(sq, id=id_, description=desc)
                seqs_for_fasta.append(srec)
    
    with open(OUT_FILE, 'w') as outf:
        # Write all the sequences as a FASTA file.
        SeqIO.write(seqs_for_fasta, outf, 'fasta')
                


    
if __name__ == "__main__":
    main(sys.argv[1:])
