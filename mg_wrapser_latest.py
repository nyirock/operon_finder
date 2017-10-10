#!/usr/bin/python

import getopt
import sys
from Bio import SeqIO
from Bio.SeqUtils import GC
import time# import time, gmtime, strftime
import os
import shutil
import pandas
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import csv
#from datetime import datetime
import numpy as np
from scipy import stats

__author__ = "Andriy Sheremet"

#Helper functions definitions

def genome_shredder(input_dct, shear_val):
    shredded = {}
    
    for key, value in input_dct.items():
        #print input_dct[i].seq
        #print i
        dic_name = key
        rec_name = value.name
        for j in range(0, len(str(value.seq)), int(shear_val)):
#            print j
            record = str(value.seq)[0+j:int(shear_val)+j]
            shredded[dic_name+"_"+str(j)] = SeqRecord(Seq(record),rec_name+"_"+str(j),'','')
            #record = SeqRecord(input_ref_records[i].seq[0+i:int(shear_val)+i],input_ref_records[i].name+"_%i"%i,"","")
    return shredded


def parse_contigs_ind(f_name):
    """
    Returns sequences index from the input files(s)
    remember to close index object after use
    """
    handle = open(f_name, "rU")
    record_dict = SeqIO.index(f_name,"fasta")
    handle.close()
    return record_dict

#returning specific sequences and overal list
def retrive_sequence(contig_lst, rec_dic):
    """
    Returns list of sequence elements from dictionary/index of SeqIO objects specific to the contig_lst parameter
    """
    contig_seqs = list()
    #record_dict = rec_dic
    #handle.close()
    for contig in contig_lst:
        contig_seqs.append(str(rec_dic[contig].seq))#fixing BiopythonDeprecationWarning
    return contig_seqs


def filter_seq_dict(key_lst, rec_dic):
    """
    Returns filtered dictionary element from rec_dic according to sequence names passed in key_lst
    """
    return { key: rec_dic[key] for key in key_lst }

def unique_scaffold_topEval(dataframe):
#returns pandas series object
    variables = list(dataframe.columns.values)
    scaffolds=dict()
    rows=list()

    for row in dataframe.itertuples():

        #if row[1]=='Ga0073928_10002560':
        if row[1] not in scaffolds:
            scaffolds[row[1]]=row
        else:
            if row[11]<scaffolds[row[1]][11]:
                scaffolds[row[1]]=row
    rows=scaffolds.values()
    #variables=['quid', 'suid', 'iden', 'alen', 'mism', 'gapo', 'qsta', 'qend', 'ssta', 'send', 'eval', 'bits']
    df = pandas.DataFrame([[getattr(i,j) for j in variables] for i in rows], columns = variables)
    return df
    
def unique_scaffold_topBits(dataframe):
#returns pandas series object
    variables = list(dataframe.columns.values)
    scaffolds=dict()
    rows=list()

    for row in dataframe.itertuples():

        #if row[1]=='Ga0073928_10002560':
        if row[1] not in scaffolds:
            scaffolds[row[1]]=row
        else:
            if row[12]>scaffolds[row[1]][12]:
                scaffolds[row[1]]=row
    rows=scaffolds.values()
    #variables=['quid', 'suid', 'iden', 'alen', 'mism', 'gapo', 'qsta', 'qend', 'ssta', 'send', 'eval', 'bits']
    df = pandas.DataFrame([[getattr(i,j) for j in variables] for i in rows], columns = variables)
    return df
    
    
def close_ind_lst(ind_lst):
    """
    Closes index objects supplied in input parameter list
    """
    for index in ind_lst:
        index.close()

def usage():
    print "\nThis is the usage function\n"
#    print 'Usage: '+sys.argv[0]+' -i <input_file> [-o <output>] [-l <minimum length>]'
#    print 'Example: '+sys.argv[0]+' -i input.fasta -o output.fasta -l 100'

def main(argv):
    
    #default parameters
    mg_lst = []
    ref_lst = []
    e_val = 1e-5
    alen = 50.0
    alen_percent = True
    alen_bp = False
    iden = 95.0
    name= "output"
    fmt_lst = ["fasta"]
    supported_formats =["fasta", "csv"]
    iterations = 1
    alen_increment = 5.0
    iden_increment = 0.0
    blast_db_Dir = ""
    results_Dir = ""
    input_files_Dir = ""
    ref_out_0 = ""
    blasted_lst = []
    continue_from_previous = False #poorly supported, just keeping the directories
    skip_blasting = False
    debugging = False
    sheared = False
    shear_val = None
    logfile = ""
    
    
             
    try:                                
        opts, args = getopt.getopt(argv, "r:m:n:e:a:i:s:f:h", ["reference=", "metagenome=", "name=", "e_value=", "alignment_length=", "identity=","shear=","format=", "iterations=", "alen_increment=", "iden_increment=","continue_from_previous","skip_blasting","debugging", "help"])
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

        elif opt in ("--continue_from_previous"):
            continue_from_previous = True
            if debugging:
                print "Continue after failure:", continue_from_previous
        elif opt in ("--debugging"):
            debugging = True
            if debugging:
                print "Debugging messages:", debugging   
                     
        elif opt in ("-r", "--reference"):
            if arg:
                ref_lst=arg.split(',')
                #infiles = arg
            if debugging:
                print "Reference file(s)", ref_lst  
        elif opt in ("-m", "--metagenome"):
            if arg:
                mg_lst=arg.split(',')
                #infiles = arg
            if debugging:
                print "Metagenome file(s)", mg_lst
            
        elif opt in ("-f", "--format"):
            if arg:
                fmt_lst=arg.split(',')
                #infiles = arg
            if debugging:
                print "Output format(s)", fmt_lst
        
        elif opt in ("-n", "--name"):
            if arg.strip():              
                name = arg
            if debugging:
                print "Project name", name 
            
        elif opt in ("-e", "--e_value"):
            try:
                e_val = float(arg)
            except:
                print "\nERROR: Please enter numerical value as -e parameter (default: 1e-5)"
                usage()
                sys.exit(1)
            if debugging:
                print "E value", e_val
            
        elif opt in ("-a", "--alignment_length"):
            if arg.strip()[-1]=="%":
                alen_bp = False
                alen_percent = True
            else:
                alen_bp = True
                alen_percent = False
                
            try:
                alen = float(arg.split("%")[0])
            except:
                print "\nERROR: Please enter a numerical value as -a parameter (default: 50.0)"
                usage()
                sys.exit(1)
            if debugging:
                print "Alignment length", alen           
            
        elif opt in ("-i", "--identity"):
            try:
                iden = float(arg)
            except:
                print "\nERROR: Please enter a numerical value as -i parameter (default: 95.0)"
                usage()
                sys.exit(1)
            if debugging:
                print "Alignment length", iden    
        elif opt in ("-s", "--shear"):
            sheared = True
            try:
                shear_val = int(arg)
            except:
                print "\nERROR: Please enter an integer value as -s parameter"
                usage()
                sys.exit(1)
            if debugging:
                print "Alignment length", iden 
        elif opt in ("--iterations"):
            try:
                iterations = int(arg)
            except:
                
                print "\nWARNING: Please enter integer value as --iterations parameter (using default: 1)"
            if debugging:
                print "Iterations: ", iterations  
            
        elif opt in ("--alen_increment"):
            
            try:
                alen_increment = float(arg)
            except:
                print "\nWARNING: Please enter numerical value as --alen_increment parameter (using default: )", alen_increment
            if debugging:
                print "Alignment length increment: ", alen_increment 
 
        elif opt in ("--iden_increment"):
            
            try:
                iden_increment = float(arg)
            except:
                print "\nWARNING: Please enter numerical value as --iden_increment parameter (using default: )", iden_increment
            if debugging:
                print "Alignment length increment: ", iden_increment 
        elif opt in ("--skip_blasting"):
            skip_blasting = True
            if debugging:
                print "Blasting step omitted; Using previous blast output."
            
    for ref_file in [x for x in ref_lst if x]:
        try:
            #
            with open(ref_file, "rU") as hand_ref:
                pass
        except:
            print "\nERROR: Reference File(s) ["+ref_file+"] doesn't exist"
            usage()
            sys.exit(1)

            
    for mg_file in [x for x in mg_lst if x]:
        try:
            #
            with open(mg_file, "rU") as hand_mg:
                pass
        except:
            print "\nERROR: Metagenome File(s) ["+mg_file+"] doesn't exist"
            usage()
            sys.exit(1) 
            
    for fmt in [x for x in fmt_lst if x]:
        if fmt not in supported_formats:
            print "\nWARNING: Output format [",fmt,"] is not supported"
            print "\tUse -h(--help) option for the list of supported formats"
            fmt_lst=["fasta"]
            print "\tUsing default output format: ", fmt_lst[0]
 
    
    project_dir = name
    if not continue_from_previous:
        if os.path.exists(project_dir):
            shutil.rmtree(project_dir)
        try:
            os.mkdir(project_dir)
        except OSError:
            print "ERROR: Cannot create project directory: " + name
            raise
    
    print "\n\t Initial Parameters:"
    print "\nProject Name: ", name,'\n'
    print "Project Directory: ", os.path.abspath(name),'\n'
    print "Reference File(s): ", ref_lst,'\n'
    if sheared:
        print "Shear Reference File(s):", str(shear_val)+"bp",'\n'
    print "Metagenome File(s): ", mg_lst,'\n'
    print "E Value: ", e_val, "\n"
    if alen_percent:
        print "Alignment Length: "+str(alen)+'%\n'
    if alen_bp:
        print "Alignment Length: "+str(alen)+'bp\n'
    print "Sequence Identity: "+str(iden)+'%\n'
    print "Output Format(s):", fmt_lst,'\n'
    if iterations > 1:
        print "Iterations: ", iterations, '\n'
        print "Alignment Length Increment: ", alen_increment, '\n'
        print "Sequence identity Increment: ", iden_increment, '\n'

    #Initializing directories
    blast_db_Dir = name+"/blast_db"
    if not continue_from_previous:
        if os.path.exists(blast_db_Dir):
            shutil.rmtree(blast_db_Dir)
        try:
            os.mkdir(blast_db_Dir)
        except OSError:
            print "ERROR: Cannot create project directory: " + blast_db_Dir
            raise

    results_Dir = name+"/results"
    if not continue_from_previous:
    
        if os.path.exists(results_Dir):
            shutil.rmtree(results_Dir)
        try:
            os.mkdir(results_Dir)
        except OSError:
            print "ERROR: Cannot create project directory: " + results_Dir
            raise

    input_files_Dir = name+"/input_files"
    if not continue_from_previous:
    
        if os.path.exists(input_files_Dir):
            shutil.rmtree(input_files_Dir)
        try:
            os.mkdir(input_files_Dir)
        except OSError:
            print "ERROR: Cannot create project directory: " + input_files_Dir
            raise

# Writing raw reference files into a specific input filename
    input_ref_records = {}
    for reference in ref_lst:
        ref_records_ind = parse_contigs_ind(reference)
        #ref_records = dict(ref_records_ind)
        input_ref_records.update(ref_records_ind)
        ref_records_ind.close()
        #input_ref_records.update(ref_records)
        
    ref_out_0 = input_files_Dir+"/reference0.fna"
    if (sheared & bool(shear_val)):
        with open(ref_out_0, "w") as handle:
            SeqIO.write(genome_shredder(input_ref_records, shear_val).values(), handle, "fasta")

            #NO NEED TO CLOSE with statement will automatically close the file
    else:
        with open(ref_out_0, "w") as handle:
            SeqIO.write(input_ref_records.values(), handle, "fasta")

# Making BLAST databases
    #output fname from before used as input for blast database creation
    input_ref_0 = ref_out_0
    title_db = name+"_db"#add iteration functionality
    outfile_db = blast_db_Dir+"/iteration"+str(iterations)+"/"+name+"_db"#change into for loop
    os.system("makeblastdb -in "+input_ref_0+" -dbtype nucl -title "+title_db+" -out "+outfile_db+" -parse_seqids")
    
# BLASTing query contigs
    if not skip_blasting:
        print "\nBLASTing query file(s):"
        for i in range(len(mg_lst)):
            
            database = outfile_db # adjust for iterations
            blasted_lst.append(results_Dir+"/recruited_mg_"+str(i)+".tab")
            start = time.time()
            os_string = 'blastn -db '+database+' -query \"'+mg_lst[i]+'\" -out '+blasted_lst[i]+" -evalue "+str(e_val)+"  -outfmt 6 -num_threads 8"
            #print os_string
            os.system(os_string)
            print "\t"+mg_lst[i]+"; Time elapsed: "+str(time.time()-start)+" seconds."
    else:
        for i in range(len(mg_lst)):
            blasted_lst.append(results_Dir+"/recruited_mg_"+str(i)+".tab")
        
        
# Parsing BLAST outputs
    blast_cols = ['quid', 'suid', 'iden', 'alen', 'mism', 'gapo', 'qsta', 'qend', 'ssta', 'send', 'eval', 'bits']
    recruited_mg=[]
    for i in range(len(mg_lst)):
        try:
            df = pandas.read_csv(blasted_lst[i] ,sep="\t", header=None)
        except:
            df = pandas.DataFrame(columns=blast_cols)
        df.columns=blast_cols
        recruited_mg.append(df)
        
#    print len(recruited_mg[0])
#    print len(recruited_mg[1])

    #creating all_records entry
#! Remember to close index objects after they are no longer needed
#! Use helper function close_ind_lst()
    all_records = []
    all_input_recs = parse_contigs_ind(ref_out_0)
    
    ##calculating GC of the reference
    if (len(all_input_recs)>1):
        ref_gc_lst = np.array([GC(x.seq) for x in all_input_recs.values()])
        ref_cnt = ref_gc_lst.size
        ref_gc_avg = np.mean(ref_gc_lst)
        ref_gc_avg_std = np.std(ref_gc_lst)
        if(len(ref_gc_lst) > 0):
            ref_gc_avg_sem = stats.sem(ref_gc_lst, axis=0)
        else:
            ref_gc_avg_sem=0
    else:
        if (debugging):
            print "Only one reference"
        ref_gc_lst = np.array([GC(x.seq) for x in all_input_recs.values()])
        ref_cnt = ref_gc_lst.size
        ref_gc_avg = np.mean(ref_gc_lst)
        ref_gc_avg_std=0
        ref_gc_avg_sem=0
    #ref_gc_avg_sem = stats.sem(ref_gc_lst, axis=0)
    
#    _ = 0
#    for key, value in all_input_recs.items():
#        _ +=1
#        if _ < 20:
#            print key, len(value)
    
    
    print "\nIndexing metagenome file(s):"
    for i in range(len(mg_lst)):
        start = time.time()
        all_records.append(parse_contigs_ind(mg_lst[i]))
        print "\t"+mg_lst[i]+" Indexed in : "+str(time.time()-start)+" seconds."

# Transforming data
    print "\nParsing recruited contigs:"
    for i in range(len(mg_lst)):
        start = time.time()
    #cutoff_contigs[dataframe]=evalue_filter(cutoff_contigs[dataframe])
        recruited_mg[i]=unique_scaffold_topBits(recruited_mg[i])
        contig_list = recruited_mg[i]['quid'].tolist()
        recruited_mg[i]['Contig_nt']=retrive_sequence(contig_list, all_records[i])
        recruited_mg[i]['Contig_size']=recruited_mg[i]['Contig_nt'].apply(lambda x: len(x))
        #recruited_mg[i]['Ref_nt']=recruited_mg[i]['suid'].apply(lambda x: all_input_recs[str(x)].seq)
        recruited_mg[i]['Ref_size']=recruited_mg[i]['suid'].apply(lambda x: len(all_input_recs[str(x)]))
        
        recruited_mg[i]['Ref_GC']=recruited_mg[i]['suid'].apply(lambda x: GC(all_input_recs[str(x)].seq))
        #recruited_mg[i]['Coverage']=recruited_mg[i]['alen'].apply(lambda x: 100.0*float(x))/min(recruited_mg[i]['Contig_size'].apply(lambda y: y),recruited_mg[i]['Ref_size'].apply(lambda z: z))
        #df.loc[:, ['B0', 'B1', 'B2']].min(axis=1)
        recruited_mg[i]['Coverage']=recruited_mg[i]['alen'].apply(lambda x: 100.0*float(x))/recruited_mg[i].loc[:,["Contig_size", "Ref_size"]].min(axis=1)
        recruited_mg[i]['Metric']=recruited_mg[i]['Coverage']*recruited_mg[i]['iden']/100.0
        try:
            recruited_mg[i]['Contig_GC']=recruited_mg[i]['Contig_nt'].apply(lambda x: GC(x))
        except:
            recruited_mg[i]['Contig_GC']=recruited_mg[i]['Contig_nt'].apply(lambda x: None)
        try:
            recruited_mg[i]['Read_RPKM']=1.0/((recruited_mg[i]['Ref_size']/1000.0)*(len(all_records[i])/1000000.0))
        except:
            recruited_mg[i]['Read_RPKM']=np.nan
        
        #recruited_mg[i] = recruited_mg[i][['quid', 'suid', 'iden', 'alen','Coverage','Metric', 'mism', 'gapo', 'qsta', 'qend', 'ssta', 'send', 'eval', 'bits','Ref_size','Ref_GC','Ref_nt','Contig_size','Contig_GC','Contig_nt']]
        recruited_mg[i] = recruited_mg[i][['quid', 'suid', 'iden', 'alen','Coverage','Metric', 'mism', 'gapo', 'qsta', 'qend', 'ssta', 'send', 'eval', 'bits','Ref_size','Ref_GC','Contig_size','Contig_GC','Read_RPKM','Contig_nt']]
        print "\tContigs from "+mg_lst[i]+" parsed in : "+str(time.time()-start)+" seconds."
   
# Here would go statistics functions and producing plots
#
#
#
#
#

# Quality filtering before outputting
    if alen_percent:
        for i in range(len(recruited_mg)):
            recruited_mg[i]=recruited_mg[i][(recruited_mg[i]['iden']>=iden)&(recruited_mg[i]['Coverage']>=alen)&(recruited_mg[i]['eval']<=e_val)]
    if alen_bp:
        for i in range(len(recruited_mg)):
            recruited_mg[i]=recruited_mg[i][(recruited_mg[i]['iden']>=iden)&(recruited_mg[i]['alen']>=alen)&(recruited_mg[i]['eval']<=e_val)]
                

#    print  len(recruited_mg[0])
#    print len(recruited_mg[1])

# Batch export to outfmt (csv and/or multiple FASTA)
    alen_str = ""
    iden_str = "_iden_"+str(iden)+"%"
    if alen_percent:
        alen_str = "_alen_"+str(alen)+"%"
    if alen_bp:
        alen_str = "_alen_"+str(alen)+"bp"

    if iterations > 1:
        prefix=name+"/results/"+name.split("/")[0]+"_iter_e_"+str(e_val)+iden_str+alen_str
    else:
        prefix=name+"/results/"+name.split("/")[0]+"_e_"+str(e_val)+iden_str+alen_str
        
    if sheared:
        prefix = prefix+'_sheared_'+str(shear_val)+"bp"
        
    prefix = prefix + "_recruited_mg_"

#initializing log file data

    logfile=name.split("/")[0]+"/results_log.csv"
    try:
        run = int(name.split("/")[-1].split("_")[-1])# using "_" less depends on the wrapper script
    except:
        if name.split("/")[-1].split("_")[-1]==name:
            run = 0
        else:
            print "Warning: Run identifier could not be written in: "+logfile
            #sys.exit(1)
            run = None
    alen_header = "Min alen"
    if alen_bp:
        alen_header = alen_header+" (bp)"
    if alen_percent:
        alen_header = alen_header+" (%)"
        
    shear_header = "Reference Shear (bp)"
    shear_log_value = 0
    if sheared:
        shear_log_value = str(shear_val)
        
    
    print "\nWriting files:"

    for i in range(len(mg_lst)):
        records= []
        if "csv" in fmt_lst:
            outfile1 = prefix+str(i)+".csv"
            recruited_mg[i].to_csv(outfile1, sep='\t')
            print str(len(recruited_mg[i]))+" sequences written to "+outfile1
        if "fasta" in fmt_lst:
            ids = recruited_mg[i]['quid'].tolist()
            #if len(ids)==len(sequences):
            for j in range(len(ids)):
                records.append(all_records[i][ids[j]])
            outfile2 = prefix+str(i)+".fasta" 
            with open(outfile2, "w") as output_handle:
                SeqIO.write(records, output_handle, "fasta")
            print str(len(ids))+" sequences written to "+outfile2
            
#Writing logfile
        
        try:
            time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        except:
            print "Warning: Time identifier could not be written in: "+logfile
        metagenome = mg_lst[i]
        #contig info
        
        rpkm_lst = np.array(recruited_mg[i]['Read_RPKM'].tolist())
        if(len(rpkm_lst) > 0):
            rpkm = np.sum(rpkm_lst)
            rpkm_std= np.std(rpkm_lst)
            rpkm_sem = np.std(rpkm_lst)*np.sqrt(len(rpkm_lst))

        else:
            rpkm = 0
            rpkm_std= 0
            rpkm_sem=0
        
        
        sizes_lst = np.array(recruited_mg[i]['Contig_size'].tolist())
        if(len(sizes_lst) > 0):
            sizes_avg = np.mean(sizes_lst)
            sizes_avg_std= np.std(sizes_lst)
            if(len(sizes_lst) > 1):
                sizes_avg_sem = stats.sem(sizes_lst, axis=0)
            else:
                sizes_avg_sem = 0
        else:
            sizes_avg = 0
            sizes_avg_std= 0
            sizes_avg_sem=0
        #sizes_avg_sem = stats.sem(sizes_lst, axis=0)
        
        alen_lst = np.array(recruited_mg[i]['alen'].tolist())
        if(len(alen_lst) > 0):
            alen_avg = np.mean(alen_lst)
            alen_avg_std = np.std(alen_lst)
            if(len(alen_lst) > 1):
                alen_avg_sem = stats.sem(alen_lst, axis=0)
            else:
                alen_avg_sem = 0
        else:
            alen_avg = 0
            alen_avg_std = 0
            alen_avg_sem=0
        #alen_avg_sem = stats.sem(alen_lst, axis=0)
        
        iden_lst = np.array(recruited_mg[i]['iden'].tolist())
        if(len(iden_lst) > 0):
            iden_avg = np.mean(iden_lst)
            iden_avg_std = np.std(iden_lst)
            if(len(iden_lst) > 1):
                iden_avg_sem = stats.sem(iden_lst, axis=0)
            else:
                iden_avg_sem = 0
        else:
            iden_avg = 0
            iden_avg_std = 0
            iden_avg_sem=0
        #iden_avg_sem = stats.sem(iden_lst, axis=0)

        gc_lst = np.array(recruited_mg[i]['Contig_GC'].tolist())
        if(len(gc_lst) > 0):
            gc_avg = np.mean(gc_lst)
            gc_avg_std = np.std(gc_lst)
            if(len(gc_lst) > 1):
                gc_avg_sem = stats.sem(gc_lst, axis=0)
            else:
                gc_avg_sem = 0
        else:
            gc_avg = 0
            gc_avg_std = 0
            gc_avg_sem=0
        
        if ref_cnt > 0:
            recr_percent = float(len(ids))/float(len(all_records[i]))*100
        else:
            recr_percent = 0.0


        
        #log_header = ['Run','Project Name','Created', 'Reference(s)','Metagenome', 'No. Contigs','No. References', alen_header, "Min iden (%)", shear_header, "Mean Contig Size (bp)","STD Contig Size", "SEM Contig Size", "Mean Contig alen (bp)","STD Contig alen", "SEM Contig alen", "Mean Contig iden (bp)","STD Contig iden", "SEM Contig iden", "Mean Contig GC (%)","STD Contig GC","SEM Contig GC","Mean Reference GC (%)","STD Reference GC","SEM Reference GC"]
        log_header = ['Run','Project Name','Created', 'Reference(s)', shear_header,'No. Ref. Sequences','Metagenome','No. Metagenome Contigs' , alen_header, "Min iden (%)",'No. Recruited Contigs','% Recruited Contigs', 'Total RPKM', 'RPKM STD', 'RPKM SEM', "Mean Rec. Contig Size (bp)","STD Rec. Contig Size", "SEM Rec. Contig Size", "Mean alen (bp)","STD alen", "SEM alen", "Mean Rec. Contig iden (bp)","STD Rec. Contig iden", "SEM Rec. Contig iden", "Mean Rec. Contigs GC (%)","STD Rec. Contig GC","SEM Rec. Contig GC","Mean Total Reference(s) GC (%)","STD Total Reference(s) GC","SEM Total Reference(s) GC"]
        #log_row = [run,name.split("/")[0],time_str, ";".join(ref_lst), metagenome, len(ids),ref_cnt, alen, iden, shear_log_value, sizes_avg,sizes_avg_std, sizes_avg_sem, alen_avg,alen_avg_std, alen_avg_sem, iden_avg,iden_avg_std, iden_avg_sem, gc_avg,gc_avg_std, gc_avg_sem,ref_gc_avg,ref_gc_avg_std, ref_gc_avg_sem]
        log_row = [run,name.split("/")[0],time_str, ";".join(ref_lst), shear_log_value,ref_cnt, metagenome,len(all_records[i]) , alen, iden,len(ids),recr_percent,rpkm, rpkm_std, rpkm_sem, sizes_avg,sizes_avg_std, sizes_avg_sem, alen_avg,alen_avg_std, alen_avg_sem, iden_avg,iden_avg_std, iden_avg_sem, gc_avg,gc_avg_std, gc_avg_sem,ref_gc_avg,ref_gc_avg_std, ref_gc_avg_sem]    
        if os.path.isfile(logfile):#file exists - appending
            with open(logfile, "a") as log_handle:
                log_writer = csv.writer(log_handle, delimiter='\t')
                log_writer.writerow(log_row)
        else:#no file exists - writing
            with open(logfile,"w") as log_handle:
                log_writer = csv.writer(log_handle, delimiter='\t')
                log_writer.writerow(log_header)
                log_writer.writerow(log_row)
            
            
    close_ind_lst(all_records)
    close_ind_lst([all_input_recs])
    
    

    #run = 0

        


        

    #all_records[i].close()# keep open if multiple iterations

#recruited_mg_1 = pandas.read_csv(out_name1 ,sep="\t", header=None)
#recruited_mg_1.columns=['quid', 'suid', 'iden', 'alen', 'mism', 'gapo', 'qsta', 'qend', 'ssta', 'send', 'eval', 'bits']

#recruited_mg_2 = pandas.read_csv(out_name2 ,sep="\t", header=None)
#recruited_mg_2.columns=['quid', 'suid', 'iden', 'alen', 'mism', 'gapo', 'qsta', 'qend', 'ssta', 'send', 'eval', 'bits']

#recruited_mg = [recruited_mg_1, recruited_mg_2]




#    blast_db_Dir = ""
#    results_Dir = ""
#    input_files_Dir = ""
    
#    parsed = SeqIO.parse(handle, "fasta")
#
#    records = list()
#
#
#    total = 0
#    processed = 0
#    for record in parsed:
#        total += 1
#        #print(record.id), len(record.seq)
#        if len(record.seq) >= length:
#            processed += 1
#            records.append(record)
#    handle.close()   
#    
#    print "%d sequences found"%(total)
#    
#    try:
#        output_handle = open(outfile, "w")
#        SeqIO.write(records, output_handle, "fasta")
#        output_handle.close()
#        print "%d sequences written"%(processed)
#    except:
#        print "ERROR: Illegal output filename"
#        sys.exit(1)
    
    
    
if __name__ == "__main__":
    main(sys.argv[1:]) 
