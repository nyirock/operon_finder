#!/usr/bin/python

__author__ = "Andriy Sheremet"


import getopt
import sys
from Bio import  SeqIO
from Bio.SeqUtils import GC
import operator


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

    
    faa_filename = params['infile']+"_converted.faa"
    #print("hehe")
    #print(type(faa_filename))
    output_handle = open(faa_filename, "w")
    #print("haha")
    for seq_record in SeqIO.parse(handle, "genbank") :
        print "Dealing with GenBank record %s" % seq_record.id
        for seq_feature in seq_record.features :
            if seq_feature.type=="CDS" :
                assert len(seq_feature.qualifiers['translation'])==1
                output_handle.write(">%s from %s\n%s\n" % (seq_feature.qualifiers['locus_tag'][0],seq_record.name,seq_feature.qualifiers['translation'][0]))
                
    output_handle.close()
    handle.close()
    
if __name__ == "__main__":
    main(sys.argv[1:])
