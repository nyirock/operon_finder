#!/usr/bin/python
import tempfile
import subprocess
import shlex
import getopt
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
#from Bio.SeqUtils import GC
import time# import time, gmtime, strftime
import os
import shutil
import pandas as pd
import csv
#from datetime import datetime
import numpy as np
import re
#from scipy import stats

__author__ = "Andriy Sheremet"

#Helper functions definitions

#input_file = '76969.assembled.faa'
model='models/amoCAB'
max_distance=3
min_operon_size=2



def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch
            
            
def joinDescriptionColumns(descr_columns):
    merged=''
    for row in descr_columns.dropna():
        if (row != None):
            #print(merged)
            merged+=row+' '
    #print(descr_columns)
    return merged.strip()

def usage():
    print("\nThis is the usage function\n")
#    print 'Usage: '+sys.argv[0]+' -i <input_file> [-o <output>] [-l <minimum length>]'
#    print 'Example: '+sys.argv[0]+' -i input.fasta -o output.fasta -l 100'


    
def extractFeatures(ids, filename):
    features=[]
    pos=0
    for record in SeqIO.parse(filename,"fasta"):
        pos+=1
        if record.name in ids:
            features.append({'qid':record.id, 'position':pos, 'descr':record.description, 'seq':str(record.seq)})
    return pd.DataFrame(features)

def extractFeatures2(ids, filename):#putting in a seq object instead of string
    features=[]
    pos=0
    for record in SeqIO.parse(filename,"fasta"):
        pos+=1
        if record.name in ids:
            features.append({'qid':record.id, 'position':pos, 'descr':record.description, 'seq':record})
    return pd.DataFrame(features)

def is_operon(x):
    return x['diff1'] or x['diff2']

def rel_coordinates(x):
    global max_distance
    
    if ((x['diff1'] < max_distance) or (x['diff2']) <max_distance):
        return min(x['diff1'],x['diff2'])
    else:
        return 0

#this function combines adjacent operons(fragments) into one    input_file = '../76969.assembled.faa'

def operonCount(lst):
    opCnt=0
    state=False
    for i in xrange(len(lst)):
        if lst[i]:
            newState=True
            if (state==False) and (newState==True):
                opCnt+=1
            lst[i]=opCnt
            state=newState
        else:
            state=False
    return lst

#this function splits adjacent operons, based on distance
#for the operons of large size it will artificially break it into pieces

def operonCount2(lst, pos, max_distance, min_operon_size):
    opCnt=0
    state=False
    position=0
    
    #catching range/xrange error in python3
    
    try:
        zrange = xrange
    except NameError:
        zrange = range
    
    for i in zrange(len(lst)):
        if lst[i]:
            newState=True
            #position=pos[i]
            if (state==False) and (newState==True):
                position=pos[i]
                opCnt+=1
            if (pos[i]>position+max_distance+min_operon_size+5):#!!! needs a fix so large operons won't be affected
                opCnt+=1
    
            lst[i]=opCnt
            state=newState
        else:
            state=False
    return lst    
    
    
    
    
    
def analyzeStructure(name_string, filter_operons={'CAB':0,'ABC':0}):#fix this for the extended version
    for key in filter_operons.keys():
    
        if key in name_string:
            filter_operons[key]=1
            return (key, 1)
        if key[::-1] in name_string:
            filter_operons[key]=-1
            return (key, -1)
    return (None, None)

def concatSeparatedValues(row, sep=" "):
    concatVal=''
    for i in range(output['type'].iloc[1].size):
        concatVal+=str(row.iloc[i])+sep
    concatVal=concatVal[:-len(sep)]#removes last unnesessary separator
    return concatVal.rstrip()

def operonParser(op_lst, input_file, op_type=False, full_operons=False):
    
    
    output_lst=list()#list of seq elements for fasta exporting
    
    if op_lst == []:
        return 1
    if full_operons:
        out= pd.DataFrame(columns=['operon #','input_file','type',  'rel_distance','size', 'concat_ids','concat_eval', 'concat_descr', 'C', 'A', 'B', 'full_operon'])
    out= pd.DataFrame(columns=['operon #','input_file','type','rel_distance','concat_eval','concat_ids','size',  'concat_descr', 'C', 'A', 'B'])
    
    for operon in op_lst:
        
        #print(operon)
        operon_structure=''.join(operon['type'].tolist())
        #print(operon_structure)
        op_type,forward=analyzeStructure(operon_structure)
        if op_type:
            #converting evalue to float
            operon[4]=operon[4].astype(float)#for correct sorting
            #reverse operon if needed
            #print(operon)
            if forward == -1:
                #''.join(operons[1].reset_index(drop=True)['type'].sort_index(ascending=False).tolist())
                operon=operon.reset_index(drop=True).sort_index(ascending=False)#rookie way   
            else:
                operon=operon.reset_index(drop=True) #else clause is a bit unnecessary
            #print(type(operon))
            cnt=int(operon['operon_count'].mean())#!!! needs an assert for integer value
            inp_file=input_file
            
            #ot=operon_structure #inclues extended understanding of an operon
            ot=op_type
            
            
            rd=''.join(operon.loc[:, 'rel_coordinates'].astype(str).tolist())
            sz=len(operon.index)
            
            
            idx_a=operon[operon.loc[:,'type']=='A'].sort_values(4, ascending=True).head(1).index[0]
            idx_b=operon[operon.loc[:,'type']=='B'].sort_values(4, ascending=True).head(1).index[0]
            idx_c=operon[operon.loc[:,'type']=='C'].sort_values(4, ascending=True).head(1).index[0]
            
            
            
            concat_eval= ";".join([str(x) for x in operon.loc[[idx_c, idx_a, idx_b]][4].astype(str)])
            #peron
            
            
            concat_ids=";".join([x.id for x in operon.loc[[idx_c, idx_a, idx_b]]['seq']])
            concat_descr=";".join([x.description for x in operon.loc[[idx_c, idx_a, idx_b]]['seq']])
            #str(operons[1][operons[1]['type']=='C']['seq'].values[0].seq)
            #takes one value from operon with the smallest e-vaue if multiple values are present
            #str(record.seq)
#             c_seq=str(operon[operon.loc[:,'type']=='C'].sort_values(4, ascending=True).head(1)['seq'].values[0].seq)
#             a_seq=str(operon[operon.loc[:,'type']=='A'].sort_values(4, ascending=True).head(1)['seq'].values[0].seq)
#             b_seq=str(operon[operon.loc[:,'type']=='B'].sort_values(4, ascending=True).head(1)['seq'].values[0].seq)
            c_seq=str(operon.loc[idx_c]['seq'].seq)
            a_seq=str(operon.loc[idx_a]['seq'].seq)
            b_seq=str(operon.loc[idx_b]['seq'].seq)
    
            record=SeqRecord(Seq(c_seq+a_seq+b_seq, IUPAC.protein), id=concat_ids, description=concat_descr)
#           records.append(record)
            
            output_lst.append(record)
#             if full_operons:
#                 full_op=
            
            
            out=out.append(pd.Series({'operon #':cnt, 
                                  'input_file':inp_file, 
                                  'type':ot, 
                                  'rel_distance':rd, 
                                  'size': sz,
                                  'concat_eval': concat_eval,
                                  'concat_ids': concat_ids,
                                  'concat_descr': concat_descr,
                                  'C':c_seq,
                                  'A':a_seq,
                                  'B':b_seq}), ignore_index=True)
            
            
    return out, output_lst  
    
    
def main(argv):
    
    #default constants
    input_file=''
    global model, max_distance, min_operon_size
    try:                                
        opts, args = getopt.getopt(argv, "i:m:h", ["input=","model=", "help"])
    except getopt.GetoptError:          
        usage()                         
        sys.exit(2)                     
    for opt, arg in opts:                
        if opt in ("-h", "--help"):      
            usage()              
            sys.exit() 
#        elif opt in ("--recover_after_failure"):
#            recover_after_failure = True
#            print "Recover after failure:", recover_after_failure  

        
                     
        elif opt in ("-i", "--input"):
            if arg:
                input_file=arg.strip()
                #infiles = arg
 
        
#            
#    
    try:
        #
        with open(input_file, "rU") as hand_ref:
            pass
    except:
        print("\nERROR: Input File ["+input_file+"] doesn't exist")
        usage()
        sys.exit(1)

            

    
    files=list()
    sizes={}
    #input_fasta="../76969.assembled.faa"
    #input_fasta="../2517287028.genes.faa"#m.rosea
    try:
        record_iter = SeqIO.parse(open(input_file),"fasta")
    except:
        print("\nERROR: Could not parse Input File ["+input_file+"]")
        
    basename=os.path.basename(input_file)
    bn_lst=basename.split(".")        


    for i, batch in enumerate(batch_iterator(record_iter, 30000)):
        
        label = i
        #f = tempfile.NamedTemporaryFile(delete=False)#exists on closing
        f = tempfile.NamedTemporaryFile(mode='w+t')#deleted after f.close()
        files.append((label,f))
        #f.seek(0)
        count = SeqIO.write(batch, f, "fasta")
        f.flush() #this solves the EOF problem
    #     with open(filename, "w") as handle:
    #         count = SeqIO.write(batch, handle, "fasta")
        #print("Wrote %i records to %s" % (count, f.name))
        sizes[f.name]=count
    
    
    
    
    
    
    #creating dictionary containing labels # may be entirely useles/redundant
    files2={j.name:i for i,j in files}
    
    
    #main heavylifting part

    #hmmscan --tblout amocab_hmmscan_mros.tab  amoCAB 2517287028.genes.faa > /dev/null
    #os.system("makeblastdb -in "+input_ref_0+" -dbtype nucl -title "+title_db+" -out "+outfile_db+" -parse_seqids")
    fnames=""
    flabels=""
    for (l,f) in files:
        fnames+=f.name+" "
        flabels+=str(l)+" "
    #print flabels
        #print(l)
        #print(f.name)
        #parallel -j $1 blastall "-p blastn -d nt/nt -i "{}" -o "{.}".xml -e 1e-10 -m 7 -K 100 -b 200" ::: *.fa
        #print('Exists after close:', os.path.exists(f.name))
    #this prints nicely to stdout, and is caught by p.stdout.read()
    #cmd="parallel -j 8 hmmscan "+"-o /dev/null --noali --cpu 1"+" --tblout >(tee /dev/stdout)  " + "../amoCAB {} ::: "+fnames
    #this prints to the file and stdout
    #0.tab should be deleted if exists
    if os.path.exists("0.tab"):
        os.remove("0.tab")
    cmd="parallel -j 8 hmmscan "+"-o /dev/null --noali --cpu 1"+" --tblout >(tee -a 0.tab)  " + model+" {} ::: "+fnames
    args = shlex.split(cmd)
    #os.system(cmd)
    p=subprocess.Popen(args,stdout=subprocess.PIPE)
    #print(p.stdout.read())
    output=p.stdout.read()
    #os.system(cmd)
    #command="parallel -j 8 hmmscan "+"--cpu 1"+" --tblout {}.tab " + "../amoCAB {} ::: "+fnames
    #print(cmd)
    
    
    #splitting output with reges, since hmmscan sometimes truncates spaces
    

    chunks=re.split(r"#\s+--- full sequence ---- --- best 1 domain ---- --- domain number estimation ----", output.decode())
    
    #sorting here might be useles
    total=[]
    df_scramb=pd.DataFrame()
    test_cnt=0
    for chunk in chunks:
        if chunk != '':
            #print((chunk))
    
            data=chunk.split('\n')[3:-11]
    #         if test_cnt <2:
    #             print(data)
            test_cnt+=1
            #print(data)
            if data:
                d=pd.DataFrame([x.split() for x in data])
            else:
                continue
            
            query=re.findall('# Query file:\s*(.*)', chunk)[0]
            label=files2[query]
            #print(query)
            #d['query file']=query
            #df.insert(loc=idx, column='A', value=new_col)
            d.insert(loc=0, column='query file', value=query)
            d.insert(loc=1, column='label', value=label)
            #d['label']=files2[query]
            df_scramb=df_scramb.append(d)
            #total+=chunk.split("\n")[3:-11]
            
    #lots of things here are just for compatibility with no-split version, 
    #they need to be replaced!!!
    #check for an empty dataframe
    if df_scramb.empty:
        print("\nNo hits found in the Input File ["+input_file+"]")
        sys.exit()
    df_scramb.sort_values(by='query file', inplace=True)#to fix performance warning
    df_scramb.set_index(['query file', 'label'], inplace=True)
    
    #df_scramb[18]
    
    df_scramb.iloc[:,18]=df_scramb.iloc[:,18:].apply(joinDescriptionColumns, axis=1)#this part is also redundant
    df_scramb=df_scramb.iloc[:,:19]
    
    
    df_scramb.sort_index(level=1, inplace=True)#sort by label and then reset index
    df_scramb.reset_index(inplace=True)
    
    
    ##hmmscan output df is combined with seq. features extracted from input_file
    #entire SeqIO.seq object is placed into 'seq' column
    names=df_scramb[2].tolist()
    df_scramb=df_scramb.merge(extractFeatures2(names, input_file), left_on=2, right_on='qid', suffixes=('',''))
    df_scramb.sort_values('position', inplace=True)
    #df_scramb.head()
    

    
    #df_scramb.sort_values('position')
    model2type={'AmoC': 'C', 'AMO': 'A', 'Monooxygenase_B': 'B'}
    
    
    
    l=df_scramb['position'].tolist()
    
    df_scramb['diff1']=abs((np.array(l)-np.array([0] + l[:-1])))#shift to left
    df_scramb['diff2']=abs(np.array(l)-np.array(l[1:]+[0]))#shift to right
    #df[(df['diff1'] <= 2) | (df['diff2'] <= 2)]
    
    df_scramb['is_operon']=df_scramb[['diff1', 'diff2']].apply(lambda x: x < max_distance).apply(is_operon, axis=1)
    
    df_scramb['rel_coordinates']=df_scramb[['diff1', 'diff2']].apply(rel_coordinates, axis=1)
    
    #df_scramb['operon_count']=operonCount(df_scramb['is_operon'].tolist(), df_scramb['position'].tolist())
    df_scramb['operon_count']=operonCount2(df_scramb['is_operon'].tolist(), df_scramb['position'].tolist(), max_distance, min_operon_size)
    
    df_scramb['type']=df_scramb[0].apply(lambda x: model2type[x])
    df_scramb.head()
    
    
    #groupping by combined operon df into operons by count and placing them into a list for processing
    operons=list()
    output=pd.DataFrame()
    for count,frame in df_scramb.groupby('operon_count', sort=False):
        if count > 0:

            operons.append(frame.copy())

    #ops=list(operons)
    
    #print(frame)
    
    if operons != []:
        final_frame, output_list=operonParser(operons, basename)
    else:
        print("\nNo operons detected in the Input File ["+input_file+"]")
        sys.exit()
    
    
    #check if dot is in the filename
    #removing the last extension
    if len(bn_lst)>1:
        outfile1="cab_"+'.'.join(bn_lst[:-1])+".tsv"
        outfile2="cab_"+'.'.join(bn_lst[:-1])+".fasta"
    else:
        outfile1="cab_"+'.'.join(bn_lst)+".tsv"
        outfile2="cab_"+'.'.join(bn_lst)+".fasta"        
    
    #outfile1 = ".tsv"
    #outfile2 = "out4.fasta"
    final_frame.to_csv(outfile1, sep='\t', header=True)
    
    #if len(ids)==len(sequences):
    records=list()
    
    with open(outfile2, "w") as output_handle:
        SeqIO.write(output_list, output_handle, "fasta")
        
        
        
        
    ####closing temporary files#####            
    for (l,f) in files:
        #print(l)
        #print(f.name)
        f.close()
                                                    
if __name__ == "__main__":
    main(sys.argv[1:]) 
