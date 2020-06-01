#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import sys
import os
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import variation
from numpy import median
import operator
from collections import Counter
from itertools import groupby
import multiprocessing
import argparse
import itertools
from difflib import SequenceMatcher
from operator import itemgetter
import gzip
pd.options.mode.chained_assignment = None

#agrparse section

docstring= """

USAGE:  python adapt_find.py <argument> 

Example: python adapt_find.py ILLUMINA

Arguments:
Valid required arguments:
(1) ILLUMINA #(for Illumina sequencing)
(2) ION_TORRENT #(for Ion Torrent sequencing)
(3) 454 #(for 454 sequencing)
(4) SOLID #(for SOLID sequencing)
Use --help for more info

DESCRIPTION
Identifies adapter sequences from raw sequencing dataset, trims, maps to the genome, if an genome index is provided

"""

parser = argparse.ArgumentParser(usage = "\n" +"python %(prog)s sequencing_platform [--output_path path/to/folder] [--input_path path/to/folder] [--files list of files] \n" +"\n" + "Usage examples:\nFor Illumina as sequencing_platform and if the current directory has all FASTQ files\npython ADAPT_find.py ILLUMINA\nFor Illumina as sequencing_platform and to specify filenames explicitly\npython ADAPT_find.py ILLUMINA --files filename1.fastq filename2.fastq\nFor Illumina as sequencing_platform and to specify path to folder containing FASTQ files\npython ADAPT_find.py ILLUMINA --input_path path/to/folder\n" + "\n" + "Description:\nIdentifies adapter sequences from raw sequencing dataset, trims, maps to the genome, if an genome index is provided\n", add_help=False)
required = parser.add_argument_group('required arguments')
required.add_argument("sequencing_platform", help = """Valid arguments:
(1) ILLUMINA #(for Illumina sequencing)
(2) SOLID #(for SOLID sequencing)
(3) ION_TORRENT #(for Ion Torrent sequencing)
(4) 454 #(for 454 sequencing)""")
optional = parser.add_argument_group('optional arguments')
optional.add_argument('--min_len', default = 15, help= 'minimum length parameter for CUTADAPT')
optional.add_argument('--max_len', default = 50, help= 'maximum length parameter for CUTADAPT')
optional.add_argument('--index',  help= 'paste path to bowtie index')
optional.add_argument('--input_path', help= 'paste path to FASTQ files')
optional.add_argument('--output_path', default = os.getcwd(), help= 'paste path to store output files')
optional.add_argument('--files', nargs='*', help= 'enter FASTQ files seperated by space')

optional.add_argument("-h", "--help", action='help', help='print help message')
args = parser.parse_args()

# get the current working directory
if (args.input_path!=None) and (args.files!=None):
   sys.exit('\nERROR: input path and files option cannnot be specified together. Only one of the two options can be specified. use ADAPT_find.py --help for more info\n%s'%(docstring))
elif (args.files!=None):
   files= args.files
   files = [f for f in files if f.split("/")[-1].endswith(".fastq") or f.split("/")[-1].endswith(".fastq.gz")]
else:
   if (args.input_path==None):
    cwd = os.getcwd()
   else:
    cwd= args.input_path
   files = [f for f in os.listdir(cwd) if f.endswith(".fastq") or f.endswith(".fastq.gz")]

if len(files)==0: 
   sys.exit('\nERROR: Could not find any FASTQ files. Please check if the input path is specified correctly. use ADAPT_find.py --help for more info\n%s'%(docstring))   

if (args.sequencing_platform!="ILLUMINA") and (args.sequencing_platform!="SOLID") and (args.sequencing_platform!="ION_TORRENT") and (args.sequencing_platform!="454"):
   sys.exit('\nERROR: Invalid argument for sequencing_platform . use ADAPT_find.py --help for more info\n%s'%(docstring))

if not os.path.exists("good-mapping"):
    os.makedirs("good-mapping")


if not os.path.exists("bad-mapping"):
    os.makedirs("bad-mapping")

    
if not os.path.exists("no_overepresented_sequences"):
    os.makedirs("no_overepresented_sequences")

if (args.index!=None):
   index = args.index


#check if blast exists (version >= 2.7)
os.system("blastn -version > vers.txt")
infile= open("vers.txt", "r")
lines = infile.readlines()
if lines != []:
   ver = float(lines[0].strip().split("blastn: ")[1][:3])
else:
   ver = "no"
if (ver == "no"):
 sys.exit('\nERROR: blast not found. Please install blast 2.7 or any version released after 2.7\n')
elif (ver < 2.7):
  sys.exit('\nERROR: adapt_find requires blast version 2.7 or higher\n')
else:
 blastn = "blastn"

#check if bowtie version exists (version >= 1.1)
if (args.index!=None):
   os.system("bowtie --version > vers.txt")
   infile= open("vers.txt", "r")
   lines = infile.readlines()
   if lines != []:
      ver1 = float(lines[0].strip().split("bowtie version")[1].split(" ")[1][:3])
   else:
      ver1 = "no"
   if (ver1 == "no"):
      sys.exit('\nERROR: bowtie not found. Please install bowtie version 1.1 or any version released after 1.1')
   elif (ver1 < 1.1):
      sys.exit('\nERROR: adapt_find requires bowtie 1.1 or higher\n')
   else:
      bowtie = "bowtie"
   
if os.popen("cutadapt --version").read().strip() == "":
   os.system('\nERROR: Cutadapt not found. Please install cutadapt')

os.system("rm vers.txt")

def worker1(f):
            filename = f.split("/")[-1].split(".")[0]
            print("processing file " + filename + ".fastq") 
            command = "mkdir "+ "-p "+ "aux_files/" + filename
            os.system(command)
            newpath = "aux_files/" + filename + "/"
            fastq_filename = filename + "_trimmed.fastq"
            query_file = newpath + filename + "_query.fa"
            subject_file = newpath + filename + "_subject.fa"
            blast_file = newpath + filename + "_blast.csv"
            blast2_file = newpath + filename + "_blast2.csv"
            filename4 = newpath + filename + "_adapters.csv"
            log = newpath + filename+ "_bowtie.txt"
            log2 = newpath + filename+ "_cutadapt.txt"
            if f.endswith(".fastq.gz"):
               fastq_lst=gzip.open(f,'rb').readlines()[1::4]
               fastq_lst = [line.strip().decode() for line in fastq_lst]
            else:
               infile= open(f)
               fastq_lst = infile.readlines()[1::4]
               fastq_lst = [line.strip() for line in fastq_lst]
            abund = len(fastq_lst)
            collapsed2 = Counter(fastq_lst)
            collapsed = sorted(collapsed2.items(), key=operator.itemgetter(1), reverse=True)
            ofile = open(query_file, "w")
            q_fil = collapsed[0][1]
            if ((abund > 5000000) and  q_fil < 14000) or ((abund > 10000000) and len(collapsed[0][0]) >76) :
               q_counter = 300
            #elif abund > 9000000:
               #q_counter = 50
            else:
               q_counter = 0
            ind_counter = 0
            for x,y in collapsed[q_counter:]:
              
              if ind_counter < 51:
                 
                 count_obj = Counter({'A':0,'T':0,'G':0,'C':0})
                 count_obj.update(x)
                 if ((min(count_obj.items(), key=itemgetter(1))[1])/len(x)) > 0.05:
                    ind_counter = ind_counter + 1
                    q_counter = q_counter + 1
                    ofile.write(">" +  str(y) +"\n" + str(x) + "\n")
                 else:
                    continue
              else:
                 break
            ofile.close()
            ofile = open(subject_file, "w")
            if ((abund > 10000000) and len(collapsed[0][0]) >76) and q_counter < 400:
               q_counter = 400
            elif (abund > 10000000) and q_counter < 100:
               q_counter = 101
            else:
               q_counter = q_counter
            ind_counter = 0 
            
            for x,y in collapsed[q_counter:]:
              if ind_counter < 201:
               
               count_obj = Counter({'A':0,'T':0,'G':0,'C':0})
               count_obj.update(x)
               if ((min(count_obj.items(), key=itemgetter(1))[1])/len(x)) > 0.05:
                    ind_counter = ind_counter + 1
                    ofile.write(">" +  str(y) +"\n" + str(x) + "\n")
               else:
                    continue
              else:
               break
            ofile.close()
            new_df = pd.DataFrame(collapsed[:100], columns=['Sequence','Count'])
            new_df2 = pd.DataFrame(collapsed[:40], columns=['Sequence','Count'])
                        
            if len(new_df) > 2:
             new_df["Sequence_length"] = [len(word) for word in new_df['Sequence']]
             new_df2["Sequence_length"] = [len(word) for word in new_df2['Sequence']]
             seq_len = new_df2["Sequence_length"].tolist()
             seq_median = new_df2["Sequence_length"].median()
             
             if (len(seq_len)) > 1:
                seq_var = variation(seq_len, axis=0)
             else:
                seq_var = 0
             if ((len(seq_len)>1) and (seq_var > 0.17) and (max(seq_len) < 50)) or ((seq_median < 30) and (max(seq_len) < 50)):
                print("There is no adapter sequence and the length distribution of top 100 collapsed reads is ..........")
                print(list(set(new_df2["Sequence_length"].tolist())))
                command = "cutadapt --trim-n -q 20 -m "+ str(args.min_len)+ " -M " + str(args.max_len) + " " +  f + " 2> " + log2 + " | cutadapt -q 20 -m "+ str(args.min_len)+ " -M " + str(args.max_len)  + " - > good-mapping/" + fastq_filename + " 2>> "+ log2
                print(command)
                os.system(command)
                if (args.index!=None):
                 command = bowtie +" --best -v 1 -p 20 " + index + " -q good-mapping/" + fastq_filename + " > check.sam 2> "+ log
                 os.system(command)
                 infile= open(log, "r")
                 lines = infile.readlines()
                 a = float(lines[1].strip().split("(")[1].split("%")[0])
                 no_reads = lines[0].strip().split("processed: ")[1]
                 if (a >=50):
                    element = [filename,"no adapter","na","na",abund,"na","na", no_reads,a,"good"]
                    return element
                 else:
                    command = "mv good-mapping/"+ fastq_filename + " bad-mapping/" + fastq_filename
                    print(command)
                    os.system(command)
                    element = [filename,"no adapter","na","na",abund,"na","na", no_reads,a,"bad"]
                    return element
                else:
                    element = [filename,"no adapter","na","na",abund,"na","na", "na","na","na"]
                    return element
             sequences = new_df["Sequence"].tolist()
             kmer = 25
             while(kmer!=0):
               kmers_collapsed = sorted(Counter([x[:kmer] for x in sequences]).items(), key=operator.itemgetter(1), reverse=True)
               if (float(kmers_collapsed[0][1])/float(len(new_df))) > 0.9:
                 break
               else: 
                 kmer = kmer - 1
             forward_adapter = kmers_collapsed[0][0]
             kmer = 25
             while(kmer!=0):
               kmers = Counter([x[-kmer:] for x in sequences])
               kmers_collapsed = sorted(Counter([x[-kmer:] for x in sequences]).items(), key=operator.itemgetter(1), reverse=True)
               if (float(kmers_collapsed[0][1])/float(len(new_df))) > 0.9:
                 break
               else: 
                 kmer = kmer - 1
             reverse_adapter = kmers_collapsed[0][0]
             
                          
             if (len(forward_adapter) > 2) and (len(reverse_adapter) > 2):
                 command = "cutadapt -q 20 -m "+ str(args.min_len)+ " -M " + str(args.max_len) + " " + "-g ^" + forward_adapter +  " " +  f + " 2> " + log2 + " | cutadapt -q 20 -m "+ str(args.min_len)+ " -M " + str(args.max_len) + " -a " + reverse_adapter + " - > good-mapping/" + fastq_filename + " 2>> "+ log2
                 print(command)
                 os.system(command)
                 infile= open(log2, "r")
                 lines = infile.readlines()
                 cut_adapt = []
                 for line in lines:
                   if "Reads with adapters" in line:
                    cut_adapt.append(line.strip().replace(" ", "").split("(")[1].split(")")[0])
                 if (args.index!=None):
                  command = bowtie +" --best -v 1 -p 20 " + index + " -q good-mapping/" + fastq_filename + " > check.sam 2> "+ log
                  print(command)
                  os.system(command)
                  adapter = "5 prime = " +forward_adapter + " & 3 prime = " + reverse_adapter
                  infile= open(log, "r")
                  lines = infile.readlines()
                  a = float(lines[1].strip().split("(")[1].split("%")[0])
                  no_reads = lines[0].strip().split("processed: ")[1]
                  if (a >=50):
                    element = [filename,"5&3'-anchored",forward_adapter,reverse_adapter,abund,cut_adapt[0],cut_adapt[1],no_reads,a,"good"]
                    return element
                  else:
                    command = "mv good-mapping/"+ fastq_filename + " bad-mapping/" + fastq_filename
                    print(command)
                    os.system(command)
                    element = [filename,"5&3'-anchored",forward_adapter,reverse_adapter,abund,cut_adapt[0],cut_adapt[1],no_reads,a,"bad"]
                    return element
                 else:
                   element = [filename,"5&3'-anchored",forward_adapter,reverse_adapter,abund,cut_adapt[0],cut_adapt[1],"na","na","na"]
                   return element
             
             #check now the median length
             command = blastn + " -task blastn-short -ungapped -max_hsps 2 -query " + query_file + " -subject " +  subject_file + " -outfmt '10 qseqid sseqid qseq length evalue qstart sstart' -out " + blast_file 
             os.system(command)
             csv = pd.read_csv(blast_file)
             csv.columns = ['query', 'subject', 'aligned_seq', 'aligned_length', 'evalue', 'qstart', 'sstart']
             fo = open(blast_file,"w")
             csv.to_csv(fo, sep = ",", index=False)
             fo.close()
              
             d = csv["sstart"].median()
             e = csv["qstart"].median()
             csv.insert(3, "adapter", "")
             nom_length=csv["aligned_length"].median()
             max_length=csv["aligned_length"].max()
             print("Median length of aligned sequences for filename - " + str(filename) + " is " + str(nom_length))
             
             if max_length > 100:
                csv = csv.loc[~(csv["aligned_length"] <= 60)]
                csv = csv.loc[~((csv["qstart"] <= 20) & (csv["sstart"] <= 20))]
                csv = csv.loc[~((csv["qstart"] >= 100) & (csv["sstart"] >= 100))]

             elif (nom_length >= 20) or (max_length > 50) or ((abund > 10000000) and len(collapsed[0][0]) >76):
              if (seq_median-nom_length) <= 10:
                csv = csv.loc[~((csv["qstart"] <= 5) & (csv["sstart"] <= 5))]
                csv = csv.loc[~(csv["aligned_length"] < 10)]
              elif (seq_median-nom_length) >= 40:
                csv = csv.loc[~((csv["qstart"] <= 5) | (csv["sstart"] <= 5))]
                csv = csv.loc[~((csv["qstart"] <= 15) & (csv["sstart"] <= 15))]
                csv = csv.loc[~(csv["aligned_length"] < 8)]
              else:
                csv = csv.loc[~((csv["qstart"] <= 5) | (csv["sstart"] <= 5))]
                csv = csv.loc[~((csv["qstart"] <= 15) & (csv["sstart"] <= 15))]
                csv = csv.loc[~(csv["aligned_length"] < 15)]
                if csv["sstart"].median() >= 15:
                   csv = csv.loc[~(csv["sstart"] < 15)]
                if csv["qstart"].median() >= 15:
                   csv = csv.loc[~(csv["qstart"] < 15)]
                csv.reset_index(drop=True,inplace=True)
             elif (10 < nom_length < 20) :
                csv = csv.loc[~((csv["qstart"] <= 5) | (csv["sstart"] <= 5))]
                csv = csv.loc[~((csv["qstart"] <= 15) & (csv["sstart"] <= 15))]
                csv = csv.loc[~(csv["aligned_length"] < 10)]
                csv.reset_index(drop=True,inplace=True)
             else:
                csv = csv.loc[~((csv["qstart"] <= 15) | (csv["sstart"] <= 15))]
                csv.reset_index(drop=True,inplace=True)
             csv.sort_values('aligned_length', ascending=False, inplace=True)
             csv.reset_index(drop=True,inplace=True)
             
             

             if (nom_length <= 50) and (max_length >= 60):
              csv = csv.loc[~(csv["aligned_length"] >= 60)]
              csv.reset_index(drop=True,inplace=True)
              
              
             if (len(csv)==0):
                    element = [filename,"no adapter","na","na",abund,"na","na","na","na","na"]
                    command = "cp "+ filename + ".fastq" + " no_overepresented_sequences/"
                    print(command)
                    os.system(command)
                    return element
     
             nasa = pd.DataFrame([],columns=list(csv))
             grouped = csv.groupby('query')
             
             def trim(group):
                 group.reset_index(drop=True,inplace=True)
                 for i in range (0,len(group)):
                     b = group.at[i,'aligned_length']
                     if (b >= 21):
                        seq = group.at[i,'aligned_seq'][:20]
                        group.at[i,'adapter']= seq
                     else:
                        seq = group.at[i,'aligned_seq']
                        group.at[i,'adapter']= seq
                 m=group["qstart"].tolist()
                 # group most_common output by frequency
                 freqs = groupby(Counter(m).most_common(), lambda x:x[1])
                 # pick off the first group (highest frequency)
                 q = [val for val,count in next(iter(freqs))[1]]
                 q.append(q[0]+1)
                 q.append(q[0]+2)
                 q.append(q[0]+3)
                 group = group.loc[group["qstart"].isin(q)] 
                 group = group.loc[group["evalue"] < 0.5]
                 group.reset_index(drop=True,inplace=True)
                 return group
             if (nom_length > 10):
              for name,group in grouped:
                 result = trim(group)
                 if len(result) !=0:
                  nasa=nasa.append(result, ignore_index=True) 
             else:
                 nasa = csv.loc[(csv["aligned_length"] < 15)] 
             if (len(nasa)==0):
                    element = [filename,"no adpater","na","na", abund, "na","na","na", "na","na"]
                    command = "cp "+ filename + " no_overepresented_sequences/"
                    print(command)
                    os.system(command)
                    return element
             print("Writing BLAST output")
             nasa.insert(4, "adapter_length", "")
             nasa['adapter_length'] = [len(word) for word in nasa['adapter']]
             if (nom_length <= 10):
                nasa = nasa.sort_values('aligned_seq', ascending=True )
                nasa.reset_index(drop=True,inplace=True)
                nasa['adapter'] = [word for word in nasa['aligned_seq']]
                i=0
                v= 0
                master_df =[]
                while(i<=(len(nasa)-2)):
                 seed = nasa.at[i,'aligned_seq']
                 match = nasa.at[i+1,'aligned_seq']
                 k = len(match)-len(seed) 
                 if (seed == match[:-k]) or (seed == match):
                   j=i+1
                   while(j<=(len(nasa))-1):
                     seed = nasa.at[i,'aligned_seq']
                     match = nasa.at[j,'aligned_seq']
                     m = len(match)-len(seed)
                     if (seed == match) or (seed == match[:-m]):
                        j = j+1
                     else:
                        break
                   df = "nasa_" + str(v)
                   df_temp = nasa.iloc[i:j,0:8]
                   df_temp.reset_index(drop=True,inplace=True)
                   df =df_temp
                   master_df.append(df)
                   v=v+1
                   i = j
                 else:
                   i = i+1
                if (len(master_df) >= 2):
                   
                   master_df = sorted(master_df, key=len)
                   master_df = master_df[-2:]
                   adap_lst = []
                   for ad in master_df[::-1]:
                       adap_lst.append(ad.at[0,'aligned_seq'])
                   adapter = adap_lst[0]
                   for i in range(0,(len(adap_lst)-1)):
                    seed = adap_lst[i+1]
                    match = SequenceMatcher(None, adapter, seed).find_longest_match(0, len(adapter), 0, len(seed))
                    temp = (adapter[match.a: match.a + match.size])
                    if len(temp) >= 5:
                     adapter = temp
                   for elem in adap_lst:
                    if (adapter == elem[:len(adapter)]):
                       if len(elem) > len(adapter):
                          adapter = elem
                   fo = open(filename4,"w")
                   nasa.to_csv(fo, sep = ",", index=False)
                   fo.close()
                   k = len(adapter)-nom_length
                   if (k > 1):
                    adapter = adapter[:-int(k)]
                elif (len(master_df) == 1):
                   max_df = max(master_df, key=len)
                   adapter = max_df.at[0,'aligned_seq']
                   print(nasa.head())
                   fo = open(filename4,"w")
                   nasa.to_csv(fo, sep = ",", index=False)
                   fo.close()
                   k = len(adapter)-nom_length
                   if (k > 1):
                    adapter = adapter[:-int(k)]
                else:
                   nasa['adapter'] = [word for word in nasa['aligned_seq']]
                   fo = open(filename4,"w")
                   nasa.to_csv(fo, sep = ",", index=False)
                   fo.close()
                   print(nasa.head())
                   master=nasa["adapter"].tolist()
                   adapter = max(master, key=master.count)
                   k = len(adapter)-nom_length
                   if (k > 1):
                    adapter = adapter[:-int(k)]
                nasa['adapter_length'] = [len(word) for word in nasa['adapter']]
                print(nasa.head())
             else:
                m=nasa["adapter_length"].median()
                #nasa = nasa.loc[nasa["adapter_length"] >= m]
                nasa['adapter'] = [word[:int(m)] for word in nasa['adapter']]
                fo = open(filename4,"w")
                nasa.to_csv(fo, sep = ",", index=False)
                fo.close()
                print(nasa.head())
                m=nasa["sstart"].tolist()
                # group most_common output by frequency
                freqs = groupby(Counter(m).most_common(), lambda x:x[1])
                # pick off the first group (highest frequency)
                q = [val for val,count in next(iter(freqs))[1]][0]
                master=nasa["adapter"].tolist()
                adap_elem = Counter(master)
                collapsed = sorted(adap_elem.items(), key=operator.itemgetter(1), reverse=True)
                lim = int(0.2 * len(nasa))
                if (len(collapsed)>3) and (lim < collapsed[1][1]):
                   collapsed = [item for item in collapsed if item[1] >= lim]
                elif (len(collapsed)>=3):
                   collapsed = collapsed[:3]
                if (q_fil < 14500) or (q > 30) or (abund > 10000000):
                 fil = 0
                 for x,y in collapsed:
                    fil = fil + y
                 if ((len(collapsed) > 2) and ((float(collapsed[1][1]/fil)) >= 0.15)) or ((abund > 10000000) and len(collapsed[0][0]) >76):
                   lim = 3
                 else:
                   lim = 2 
                 adap_lst = []
                 for x,y in collapsed[:lim]:
                    adap_lst.append(x)
                 if (len(adap_lst) > 1):
                   adapter = adap_lst[0]
                   for i in range(1,len(adap_lst)):
                    seed = adap_lst[i]
                    match = SequenceMatcher(None, adapter, seed).find_longest_match(0, len(adapter), 0, len(seed))
                    temp = (adapter[match.a: match.a + match.size])
                    if len(temp) > 10:
                       adapter = temp
                   for elem in adap_lst:
                    if (adapter == elem[:len(adapter)]):
                       if len(elem) > len(adapter):
                          adapter = elem
                 else:
                  adapter = collapsed[0][0]
                else:
                   adapter = collapsed[0][0]
             
             print ("\n" + "Putative three prime end adapters for filename - " +  filename + " is ")
             print (adap_lst)
             if len(forward_adapter) > 3:
                print("Five prime end adapter sequence for filename - " +  filename + " is " + str(forward_adapter))
             print("Three prime end adapter sequence for filename - " +  filename + " is " + str(adapter) + "\n" + "Trimming with CUTADAPT" + "\n")

             if len(forward_adapter) > 3:
                 command = "cutadapt -g ^" + forward_adapter + " -m " + str(args.min_len) + " -M " + str(args.max_len) + " " +  f + " 2> " + log2 + " | cutadapt -q 20 " + "-m " + str(args.min_len)+ " -M " + str(args.max_len) + " -a " + adapter + " - > good-mapping/" + fastq_filename + " 2>> "+ log2
                 print(command)
                 os.system(command)
                 infile= open(log2, "r")
                 lines = infile.readlines()
                 cut_adapt = []
                 for line in lines:
                   if "Reads with adapters" in line:
                    cut_adapt.append(line.strip().replace(" ", "").split("(")[1].split(")")[0])
                 if (args.index!=None):
                  command = bowtie +" --best -v 1 -p 20 " + index + " -q good-mapping/" + fastq_filename + " > check.sam 2> "+ log
                  print(command)
                  os.system(command)
                  infile= open(log, "r")
                  lines = infile.readlines()
                  a = float(lines[1].strip().split("(")[1].split("%")[0])
                  no_reads = lines[0].strip().split("processed: ")[1]
                  if (a >=50):
                    element = [filename,"5'-anchored&3'-normal", forward_adapter,adapter,abund, cut_adapt[0],cut_adapt[1],no_reads,a,"good"]
                    return element
                  else:
                    command = "mv good-mapping/"+ fastq_filename + " bad-mapping/" + fastq_filename
                    element = [filename,"5'-anchored&3'-normal", forward_adapter,adapter,abund, cut_adapt[0],cut_adapt[1],no_reads,a,"bad"]
                    return element
                 else:
                    element = [filename,"5'-anchored&3'-normal", forward_adapter,adapter,abund, cut_adapt[0],cut_adapt[1],"na","na","na"]
                    return element
             else:
                 command = "cutadapt -q 20 -m "+ str(args.min_len)+ " -M " + str(args.max_len) + " -a " + adapter + " -o good-mapping/" + fastq_filename +  " " + f + " > "+ log2
                 print(command)
                 os.system(command)
                 def cut(x):
                  infile= open(x, "r")
                  lines = infile.readlines()
                  for line in lines:
                   if "Reads with adapters" in line:
                    valu = float(line.strip().replace(" ", "").split("(")[1].split(")")[0].split("%")[0])
                  return valu
                 d = cut(log2)
                 if (d < 25):
                  new_lst = []
                  collapsed = sorted(collapsed2.items(), key=operator.itemgetter(1), reverse=True)
                  new_lst = []
                  for x,y in collapsed:
                   count_obj = Counter({'A':0,'T':0,'G':0,'C':0})
                   count_obj.update(x)
                   if (((min(count_obj.items(), key=itemgetter(1))[1])/len(x)) >= 0.10) and count_obj['N']==0:
                      new_lst.append((x,y))
                  print ("len of new_lst is " +   str(filename) + " "  + str(len(new_lst)))
                  filtered = [(x,y) for x,y in new_lst if y==1]
                  print ("len of filtered list is " +   str(filename) + " "  + str(len(filtered)))
                  ofile = open(query_file, "w")
                  for x,y in filtered[:1000]:
                    ofile.write(">" +  str(y) +"\n" + str(x) + "\n")
                  ofile.close()
                  ofile = open(subject_file, "w")
                  for x,y in filtered[1000:5000]:
                   ofile.write(">" +  str(y) +"\n" + str(x) + "\n")
                  ofile.close()
                  command = "blastn -task blastn-short -ungapped -max_hsps 1 -query " + query_file + " -subject " +  subject_file + " -outfmt '10 qseqid sseqid qseq length evalue qstart sstart' -out " + blast2_file
                  os.system(command)
                  csv = pd.read_csv(blast2_file)
                  csv.columns = ['query', 'subject', 'aligned_seq', 'aligned_length', 'evalue', 'qstart', 'sstart']
                  csv.to_csv(blast2_file, sep = ",", index=False)
                  csv.insert(3,"adapter","")
                  csv['adapter'] = [word[:13] for word in csv['aligned_seq']]
                  #asan = csv.loc[csv["sstart"] >= 25]
                  #asan = asan.loc[asan["qstart"] >= 25]
                  asan = csv.loc[csv["aligned_length"] < 19]
                  asan = asan.loc[asan["aligned_length"] > 13]
                  asan = asan.loc[asan["evalue"] < 0.12]
                  asan = asan.loc[((asan["qstart"] >= 30) & (asan["sstart"] >= 30))]
                  fo = open(filename4,"w")
                  asan.to_csv(fo, sep = ",", index=False)
                  fo.close()
                  print (asan.head())
                  adap_lst = sorted(Counter(asan["adapter"].tolist()).items(), key=operator.itemgetter(1), reverse=True)
                  if len(adap_lst) > 10:
                     adap_lst = adap_lst[:10]
                  adap_lst2 = [ x for x,y in adap_lst]
                  if len(adap_lst2) > 1:
                   adap_lst = adap_lst2[:2]
                   print (adap_lst)
                   adapter = adap_lst2[0]
                   for i in range(1,len(adap_lst)):
                    seed = adap_lst[i]
                    match = SequenceMatcher(None, adapter, seed).find_longest_match(0, len(adapter), 0, len(seed))
                    temp = (adapter[match.a: match.a + match.size])
                    if len(temp) > 10:
                       adapter = temp
                   for elem in adap_lst2:
                    if (adapter == elem[:len(adapter)]):
                     if len(elem) > len(adapter):
                       adapter = elem
                   print ("adpater for filename " + str(filename) + " is " + str(adapter))
                   command = "cutadapt -q 20 -m "+ str(args.min_len)+ " -M " + str(args.max_len) + " -a " + adapter + " -o good-mapping/" + fastq_filename +  " " + f + " > "+ log2
                   print (command)
                   os.system(command)
                   d = cut(log2)
                 d = str(d) + "%"
                 if (args.index!=None):
                  command = bowtie +" --best -v 1 -p 20 " + index + " -q good-mapping/" + fastq_filename + " > check.sam 2> "+ log
                  print(command)
                  os.system(command)
                  infile= open(log, "r")
                  lines = infile.readlines()
                  a = float(lines[1].strip().split("(")[1].split("%")[0])
                  no_reads = lines[0].strip().split("processed: ")[1]
                  if (a >=50):
                    element = [filename,"3'-normal","na",adapter,abund,"na",str(d),no_reads,a,"good"]
                    return element
                  else:
                    command = "mv good-mapping/"+ fastq_filename + " bad-mapping/" + fastq_filename
                    print(command)
                    os.system(command)
                    element = [filename,"3'-normal","na",adapter,abund,"na",str(d),no_reads,a,"bad"]
                    return element
                 else:
                    element = [filename,"3'-normal","na",adapter,abund,"na",str(d),"na","na","na"]
                    return element 
def worker2(f):
            filename = f.split(".")[0]
            command = "mkdir "+ "-p "+ "aux_files/" + filename
            os.system(command)
            newpath = "aux_files/" + filename + "/"
            fastq_filename = filename + "_trimmed.fastq"
            query_file = newpath + filename + "_query.fa"
            subject_file = newpath + filename + "_subject.fa"
            over_rep = newpath + filename + "_overrep.csv"
            over_rep2 = newpath + filename + "_overrep2.csv"
            blast_file = newpath + filename + "_blast.csv"
            filename4 = newpath + filename + "_adapters.csv"
            log = newpath + filename+ "_bowtie.txt"
            log2 = newpath + filename+ "_cutadapt.txt"
            if f.endswith(".fastq.gz"):
               fastq_lst=gzip.open(f,'rb').readlines()[1::4]
               fastq_lst = [line.strip().decode() for line in fastq_lst]
            else:
               infile= open(f)
               fastq_lst = infile.readlines()[1::4]
               fastq_lst = [line.strip() for line in fastq_lst]
            abund  = len(fastq_lst)
            collapsed = Counter(fastq_lst)
            collapsed = sorted(collapsed.items(), key=operator.itemgetter(1), reverse=True)
            ofile = open(query_file, "w")
            if len(collapsed) > 365:
             q_counter = 200
            else:
             q_counter = 0 
            ind_counter = 0 
            print("Total length of unique sequences is " + str(len(collapsed)))
            for x,y in collapsed[q_counter:]:
              q_counter = q_counter + 1
              if ind_counter < 15:
                 count_obj = Counter({'A':0,'T':0,'G':0,'C':0})
                 count_obj.update(x)
                 if ((count_obj[min(count_obj)])/len(x) > 0.01) or (len(collapsed)<365):
                    ind_counter = ind_counter + 1
                    ofile.write(">" +  str(y) +"\n" + str(x) + "\n")
                 else:
                    continue
              else:
                 break
            ofile.close()
            ofile = open(subject_file, "w")
            if (len(collapsed)> 500) and (q_counter <= 300):
               q_counter = 300
            ind_counter = 0 
            for x,y in collapsed[q_counter:]:
              if ind_counter < 100:
               
               count_obj = Counter({'A':0,'T':0,'G':0,'C':0})
               count_obj.update(x)
               if ((count_obj[min(count_obj)])/len(x) > 0.01) or (len(collapsed)<365):
                    ind_counter = ind_counter + 1
                    ofile.write(">" +  str(y) +"\n" + str(x) + "\n")
               else:
                    continue
              else:
               break
            ofile.close()
            
            new_df2 = pd.DataFrame(collapsed[:15], columns=['Sequence','Count'])
            if len(new_df2) > 2:
             
             new_df2["Sequence_length"] = [len(word) for word in new_df2['Sequence']]
             seq_len = new_df2["Sequence_length"].tolist()
             seq_median = new_df2["Sequence_length"].median()

             if (len(seq_len)) > 1:
                seq_var = variation(seq_len, axis=0)
             else:
                seq_var = 0
             if ((len(seq_len)>1) and (seq_var > 0.17) and (max(seq_len) < 50)) or ((seq_median < 30) and (max(seq_len) < 50)):
                print("There is no adapter sequences and the length distribution of sequences is ..........")
                print(seq_len)
                command = "cutadapt --trim-n" + " -m "+ str(args.min_len)+ " -M " + str(args.max_len) + " " +  f + " 2> " + log2 + " | cutadapt -q 20 -m "+ str(args.min_len)+ " -M " + str(args.max_len) + " "  + " - > good-mapping/" + fastq_filename + " 2>> "+ log2
                print(command)
                os.system(command)
                if (args.index!=None):
                 command = bowtie + " --best -v 1 -p 20 " + index + " -q good-mapping/" + fastq_filename + " > check.sam 2> "+ log
                 os.system(command)
                 infile= open(log, "r")
                 lines = infile.readlines()
                 a = float(lines[1].strip().split("(")[1].split("%")[0])
                 no_reads = lines[0].strip().split("processed: ")[1]
                 if (a >=50):
                    element = [filename,"na","na","na",abund,"na","na", no_reads,a,"good"]
                    return element
                    
                 else:
                    command = "mv good-mapping/"+ fastq_filename + " bad-mapping/" + fastq_filename
                    print(command)
                    os.system(command)
                    element = [filename,"na","na","na",abund,"na","na", no_reads,a,"bad"]
                    return element
                else:
                    element = [filename,"na","na","na",abund,"na","na","na","na","na"]
                    return element        
             #check now the median length
             command = blastn + " -task blastn-short -ungapped -max_hsps 2 -query " + query_file + " -subject " +  subject_file + " -outfmt '10 qseqid sseqid qseq length evalue qstart sstart' -out " + blast_file 
             print (command)
             os.system(command)
             csv = pd.read_csv(blast_file)
             csv.columns = ['query', 'subject', 'aligned_seq', 'aligned_length', 'evalue', 'qstart', 'sstart']
             fo = open(blast_file,"w")
             csv.to_csv(fo, sep = ",", index=False)
             fo.close() 
             df1 = csv.loc[csv["qstart"] < 15]
             df2 = csv.loc[csv["qstart"] >= 15]
             master=df1["aligned_seq"].tolist()
             adap_elem = Counter(master)
             collapsed = sorted(adap_elem.items(), key=operator.itemgetter(1), reverse=True)
             collapsed = [x for x in collapsed if len(x[0])>17]
             collapsed = [x for x in collapsed if len(x[0])<20]
             adap_lst = []
             for x,y in collapsed[:3]:
                 adap_lst.append(x)
             adapter =adap_lst[0]
             for elem in adap_lst[1:]:
                    seed = elem
                    match = SequenceMatcher(None, adapter, seed).find_longest_match(0, len(adapter), 0, len(seed))
                    temp = (adapter[match.a: match.a + match.size])
                    if len(temp) > 12:
                     adapter = temp

             fiveprime_adapter = adapter
             master=df2["aligned_seq"].tolist()
             adap_elem = Counter(master)
             collapsed = sorted(adap_elem.items(), key=operator.itemgetter(1), reverse=True)
             collapsed = [x for x in collapsed if len(x[0]) > 17]
             adap_lst = []
             for x,y in collapsed[:3]:
                 adap_lst.append(x)
             adapter =adap_lst[0]
             for elem in adap_lst[1:]:
                    seed = elem
                    match = SequenceMatcher(None, adapter, seed).find_longest_match(0, len(adapter), 0, len(seed))
                    temp = (adapter[match.a: match.a + match.size])
                    if len(temp) > 12:
                     adapter = temp
             for elem in adap_lst:
                    if (adapter == elem[:len(adapter)]):
                       if len(elem) > len(adapter):
                          adapter = elem
             threeprime_adapter = adapter                          
             d = df1["qstart"].mode()[0]
             e = df2["qstart"].mode()[0]
             if d in range(0,6):
                print("5 prime adapter OK")
             else:
                print("WARNING: Unusual query start position for fiveprime_adapter")
             if e in range(0,6):
                print("WARNING: Unusual query start position for threeprime_adapter")
             else:
                print("3 prime adapter OK")
             
             print("5prime_adapter is " + fiveprime_adapter)
             print("3prime_adapter is " + threeprime_adapter)
             adapter = "5 prime = " +fiveprime_adapter + " & 3 prime = " + threeprime_adapter
  
             command = "cutadapt -q 20 -e 0.25 -m "+ str(args.min_len)+ " -M " + str(args.max_len) + " -g " + fiveprime_adapter +  " " +  f + " 2> " + log2 + " | cutadapt -q 20 -m "+ str(args.min_len)+ " -M " + str(args.max_len) + " -a " + threeprime_adapter + " - > good-mapping/" + fastq_filename + " 2>> "+ log2
             print(command)
             os.system(command)
             infile= open(log2, "r")
             lines = infile.readlines()
             cut_adapt = []
             for line in lines:
                   if "Reads with adapters" in line:
                    cut_adapt.append(line.strip().replace(" ", "").split("(")[1].split(")")[0])
             if (args.index!=None):
              command = bowtie +" --best -v 1 -p 20 " + index + " -q good-mapping/" + fastq_filename + " > check.sam 2> "+ log
              print(command)
              os.system(command)
              infile= open(log, "r")
              lines = infile.readlines()
              a = float(lines[1].strip().split("(")[1].split("%")[0])
              no_reads = lines[0].strip().split("processed: ")[1]
              if (a >=50):
                    element = [filename,"both",fiveprime_adapter,threeprime_adapter,abund,cut_adapt[0],cut_adapt[1],no_reads,a,"good"]
                    return element
              else:
                    command = "mv good-mapping/"+ fastq_filename + " bad-mapping/" + fastq_filename
                    print(command)
                    os.system(command)
                    element = [filename,"both",fiveprime_adapter,threeprime_adapter,abund,cut_adapt[0],cut_adapt[1], no_reads,a,"bad"]
                    return element
             else: 
                    element = [filename,"both",fiveprime_adapter,threeprime_adapter,abund,cut_adapt[0],cut_adapt[1],"na","na","na"]
                    return element

def worker3(f):
    if f.endswith(".fastq"):
       os.system("perl -X csfq2fq.pl " + f + " > q_fastq/" + f.split(".fastq")[0] + "_q.fastq")
    else:
       os.system("perl -X csfq2fq.pl " + f + " > q_fastq/" + f.split(".fastq.gz")[0] + "_q.fastq")
    filename = f.split(".")[0]
    print("processing file " + f)
    command = "mkdir "+ "-p "+ "aux_files/" + filename
    os.system(command)
    newpath = "aux_files/" + filename + "/"
    query_file = newpath + filename + "_query.fa"
    subject_file = newpath + filename + "_subject.fa"
    blast_file = newpath + filename + "_blast.csv"
    infile= open("q_fastq/" + f.split(".fastq")[0] + "_q.fastq")
    fastq_lst = infile.readlines()[1::4]
    fastq_lst = [line.strip() for line in fastq_lst]
    abund = len(fastq_lst)
    collapsed = Counter(fastq_lst)
    collapsed = sorted(collapsed.items(), key=operator.itemgetter(1), reverse=True)
    ofile = open(query_file, "w")
    q_counter = 0
    ind_counter = 0
    for x,y in collapsed[q_counter:]:
       if ind_counter < 51:
          count_obj = Counter({'A':0,'T':0,'G':0,'C':0})
          count_obj.update(x)
          if ((min(count_obj.items(), key=itemgetter(1))[1])/len(x)) > 0.05:
             ind_counter = ind_counter + 1
             q_counter = q_counter + 1
             ofile.write(">" +  str(y) +"\n" + str(x) + "\n")
          else:
             continue
    ofile.close()
    ind_counter = 0
    ofile = open(subject_file, "w")
    for x,y in collapsed[q_counter:]:
       if ind_counter < 201:
          count_obj = Counter({'A':0,'T':0,'G':0,'C':0})
          count_obj.update(x)
          if ((min(count_obj.items(), key=itemgetter(1))[1])/len(x)) > 0.05:
             ind_counter = ind_counter + 1
             q_counter = q_counter + 1
             ofile.write(">" +  str(y) +"\n" + str(x) + "\n")
          else:
             continue
    ofile.close()
    command = blastn + " -task blastn-short -ungapped -max_hsps 2 -query " + query_file + " -subject " +  subject_file + " -outfmt '10 qseqid sseqid qseq length evalue qstart sstart' -out " + blast_file 
    os.system(command)
    csv = pd.read_csv(blast_file)
    csv.columns = ['query', 'subject', 'aligned_seq', 'aligned_length', 'evalue', 'qstart', 'sstart']
    #csv  =csv.loc[csv["aligned_length"] >=15]
    adapter = sorted(Counter(csv["aligned_seq"].tolist()).items(), key=operator.itemgetter(1), reverse=True)[0][0][:21]
    print("Predicted adapter sequence for filename " + filename + " is " + adapter)
    print("Trimming ~" + filename)
    newpath = "aux_files/" + filename + "/"
    log = newpath + filename + "_bowtie.txt"
    log2 = newpath + filename + "_cutadapt.txt"
    tsv_file = newpath + filename+ "_sam.tsv"
    fastq_filename = "solid-adapter-trimmed/"+ filename + "_trimmed.fastq"
    sam = newpath + filename + "_mapped.sam"
    if float(os.popen("cutadapt --version").read().strip()) > 1.18:
       sys.exit('\nERROR: Cutadapt version should be less than or equal to 1.18. Colorspace reads are supported only in cutadapt version 1.18 or earlier \n')
    command = "cutadapt -c --format=sra-fastq -a " + adapter + " -q 20 -m 15 -M50 " + filename +  ".fastq 2> " + log2 + " | cutadapt -c -q 20 -m 15 -M50 - > "  +  fastq_filename + " 2>> "+ log2
    print(command)
    os.system(command)
    infile= open(log2, "r")
    lines = infile.readlines()
    cut_adapt = []
    for line in lines:
        if "Reads with adapters" in line:
           cut_adapt.append(line.strip().replace(" ", "").split("(")[1].split(")")[0])
    if (args.index!=None):
      index = args.index
      command = bowtie +" --best -C -v 3 -p 20 " + index + " -q " + fastq_filename + " -S 2> " + log + "| samtools view -Sh -F 4 - > " + sam
      print(command)
      os.system(command)
      os.system("rm " + fastq_filename)
      os.system("egrep -v '@HD|@SQ|@PG' " + sam + " > " + tsv_file)
      tsv = pd.read_csv(tsv_file,sep='\t', header = None, quoting=3)
      if len(list(tsv)) == 16:
         tsv.columns = ['header', 'flag', 'reference', 'sequence_start', 'MAPQ','CIGAR', 'REf-MATE', 'REF-MATE-POS', 'INSERT-SIZE', 'READ', 'ASCII', 'OPT-1','OPT-2','OPT-3','OPT-4','OPT-5']
      else:
         tsv.columns = ['header', 'flag', 'reference', 'sequence_start', 'MAPQ','CIGAR', 'REf-MATE', 'REF-MATE-POS', 'INSERT-SIZE', 'READ', 'ASCII', 'OPT-1','OPT-2','OPT-3','OPT-4']
      seq = []
      ofile = open(fastq_filename, "w")
      seq = tsv[['READ', 'ASCII']].apply(tuple, axis=1).tolist()
      i = 0
      for x,y in seq:
         ofile.write(">" +  filename + "-" + str(i) +"\n" + str(x) + "\n"+ "+" +"\n"+ str(y)+"\n")
         i = i + 1
      ofile.close()
      infile= open(log, "r")
      lines = infile.readlines()
      a = float(lines[1].strip().split("(")[1].split("%")[0])
      no_reads = lines[0].strip().split("processed: ")[1]
      os.system(command)
      if (a >=50):
        return [filename,"3prime","na",adapter,abund,"na",cut_adapt[0],no_reads,a,"good"]
      else:
        command = "mv good-mapping/"+ fastq_filename + " bad-mapping/" + fastq_filename
        print(command)
        os.system(command)
        return [filename,"3prime","na",adapter,abund,"na",cut_adapt[0],no_reads,a,"bad"]   
    else:
        return [filename,"3prime","na",adapter,abund,"na",cut_adapt[0],"na","na","na"]
if __name__ == "__main__":
 if (args.sequencing_platform == "ILLUMINA"):
   worker = worker1
 elif (args.sequencing_platform == "454") or (args.sequencing_platform == "ION_TORRENT"):
   worker = worker2
 elif (args.sequencing_platform == "SOLID"):
  if not os.path.exists("q_fastq"):
    os.makedirs("q_fastq")
  if not os.path.exists("solid-adapter-trimmed"):
    os.makedirs("solid-adapter-trimmed")
  if (args.index!=None):
   os.system("samtools --version > vers.txt")
   infile= open("vers.txt", "r")
   lines = infile.readlines()
   if lines != []:
    ver = lines[0].strip().split("samtools ")[1]
   else:
    ver = "no"
   if (ver=="no"):
     sys.exit('\nERROR: SAMtools not found. Please install samtools \n')
  os.system("wget https://gist.github.com/pcantalupo/9c30709fe802c96ea2b3/archive/b5a290a3993a4845d3766a018837557bd0f0047b.zip")
  os.system("unzip -j b5a290a3993a4845d3766a018837557bd0f0047b.zip 9c30709fe802c96ea2b3-b5a290a3993a4845d3766a018837557bd0f0047b/csfq2fq.pl")
  os.system("rm -r b5a290a3993a4845d3766a018837557bd0f0047b.zip")
  worker = worker3
 else:
   sys.exit('\nERROR: Invalid argument. \n%s'%(docstring))
 size = {}
 for f in files:
  size[f] = os.path.getsize(f)

 size = sorted(size.items(), key=operator.itemgetter(1))
 asan = [k for k, v in size if v >= 10000000000]
 nasa = [k for k, v in size if v < 10000000000]
 print("Total number of input files = " + str(len(files)) + "\n" + " Number of bigger files (>=10GB) = " + str(len(asan)) + "\n" + " Number of smaller files (<10 GB) = " + str(len(nasa)))

 result_list =[]

 if (len(nasa) != 0):
    cp_count = multiprocessing.cpu_count()
    print("number of CPU is " + str(cp_count))
    if ((len(nasa)/cp_count) < 2):
       pool = multiprocessing.Pool(processes = int(cp_count/2))
       print("Processing " + str(len(nasa)) + " files")
       result_list.append(pool.map(worker, [f for f in nasa ]))
       pool.close()
       pool.join()
    else:
     denom  = int(len(nasa)/(cp_count/2))
     new_a2 = int(len(nasa)/denom)
     pool = multiprocessing.Pool(processes = int(cp_count/2))
     print("Working on first batch")
     result_list.append(pool.map(worker, [f for f in nasa[:new_a2] ]))
     pool.close()
     pool.join()
     new_a = new_a2 
     for no in range(2,(denom+1)):
      pool = multiprocessing.Pool(processes = int(cp_count/2))
      print("Working on batch#" + str(no))
      result_list.append(pool.map(worker, [f for f in nasa[new_a:(new_a + new_a2)] ]))
      new_a = new_a +new_a2
      pool.close()
      pool.join()
     #pool = multiprocessing.Pool(processes = (cp_count-5))
     if len(nasa) != (new_a-1):
      pool = multiprocessing.Pool(processes = 5)
      print("Working on batch#" + str(no+1))
      result_list.append(pool.map(worker, [f for f in nasa[(new_a):] ]))
      pool.close()
      pool.join()
 if (len(asan) != 0):
    if len(asan) >= 10:
     denom = int(len(asan)/5)
     new_a2 = int(len(asan)/denom)
     print("No of big files is " + str(len(asan)))
     pool = multiprocessing.Pool(processes = 3)
     print("Working on first batch")
     result_list.append(pool.map(worker, [f for f in asan[:new_a2] ]))
     pool.close()
     pool.join()
     new_a = new_a2
     for no in range(2,(denom+1)):
      pool = multiprocessing.Pool(processes = 5)
      print("Working on batch#" + str(no))
      result_list.append(pool.map(worker, [f for f in asan[new_a:(new_a + new_a2)] ]))
      new_a = new_a +new_a2
      pool.close()
      pool.join()
     if len(asan) != (new_a-1):
      pool = multiprocessing.Pool(processes = 5)
      print("Working on batch#" + str(no+1))
      result_list.append(pool.map(worker, [f for f in asan[new_a:] ]))
      pool.close()
      pool.join()
    else:
      print("Working on bigger file(s) (file size greater than 10 GB)")
      pool = multiprocessing.Pool(processes = 5)
      result_list.append(pool.map(worker, [f for f in asan ]))
      pool.close()
      pool.join()
 result_list = list(itertools.chain.from_iterable(result_list))  
 asan  = pd.DataFrame.from_records(result_list, columns=['filename', 'adapter_type','5prime-adapter','3prime-adapter','#raw_reads','%reads_5prime-adapter','%reads_3prime-adapter','#reads_after_trimming','mapping-%','status'])       
 asan.to_csv("master_adapters.csv", sep = ",", index=False)

 if len(files) != len(asan):
   print(str(len(files)-len(asan)) + " files have skipped trimming process")

 if (args.index!=None) and (args.sequencing_platform != "SOLID"):
   command = "rm check.sam"
   #print(command)
   os.system(command)

if (args.sequencing_platform == "SOLID"):
   os.system("rm csfq2fq.pl")
