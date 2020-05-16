#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from builtins import str
import sys
import os
import pandas as pd
import multiprocessing
from collections import Counter
import itertools
import operator
import argparse

docstring= """

USAGE:
python random_mer.py

Arguments:
Valid arguments:

Requied argument:
(1) PATH to genome file

Optional argument:
(2) --mapping_stats = "YES" or "NO"
(3) --index = PATH to bowtie index (applicable only if --mapping_stats= "YES")
(4) --input_path = PATH to FASTQ files (default current working directory)
(5) --output_path = PATH to output folder to store output file(If specified PATH doesn't exist, it will be created.Default - current working directory)
(6) --files = enter FASTQ files seperated by space(if files are not in current directory,enter the absolute path)
Use --help for more info

Argument option for --mapping_stats is  case sensitive

DESCRIPTION
This tool identifis random N-mers from one or both sides of read sequences. Trims random N-mers using cutadapt

"""



parser = argparse.ArgumentParser(usage = "\n" +"python %(prog)s ref_genome [--index path/to/index][--input_path path/to/folder][--output_path path/to/folder] [--files list of files] \n" +"\n" + "\n" + "Description:\nIdentify and trim random nucleotides from sequencing reads\n", add_help=False)
required = parser.add_argument_group('required arguments')
required.add_argument("ref_genome", help = "enter the path to the FASTA file for reference genome")
optional = parser.add_argument_group('optional arguments')
optional.add_argument('--mapping_stats', default = "NO", help = 'yes or no, if trimmed reads are to be mapped to the genome')
optional.add_argument('--index', help= 'enter path to BOWTIE index')
optional.add_argument('--input_path', default = os.getcwd(), help= 'paste path to FASTQ files')
optional.add_argument('--output_path', default = os.getcwd(), help= 'paste path to store output files')
optional.add_argument('--files', nargs='*', help= 'enter FASTQ files seperated by space')


optional.add_argument("-h", "--help", action='help', help='print help message')
args = parser.parse_args()

subject_file = args.ref_genome
if (len(subject_file.split(".fa")) ==1) and (len(subject_file.split(".fasta")) == 1) and (len(subject_file.split(".FASTA")) == 1):
    sys.exit('\nERROR:Input file should be in FASTA format\n%s'%(docstring)) 

if (args.mapping_stats == "YES") and (args.index == None):
      sys.exit('\nERROR: If mapping_stats value is YES, then the index has to be specified for bowtie mapping. use random_mer.py --help for more info\n%s'%(docstring))

if (args.mapping_stats != "YES") and (args.mapping_stats != "NO"):
   sys.exit('\nERROR: Invalid argument. use random_mer.py --help for more info\n%s'%(docstring))

if args.index != None:
   index = args.index
if args.files == None:
 os.chdir(args.input_path)
 #get all fastq files
 cwd= os.getcwd()
 files = [f for f in os.listdir(cwd) if f.endswith(".fastq")]
else:
 files = args.files

if (len(files) == 0):
  sys.exit('\nERROR: No fastq files found  \n%s' %(docstring))

os.chdir(args.output_path)


#check if Package: blast 2.7.1 exists
os.system("blastn -version > vers.txt")
infile= open("vers.txt", "r")
lines = infile.readlines()
if lines != []:
   ver = lines[0].strip().split("blastn: ")[1]
else:
   ver = "no"
if (ver != "2.7.1+"):
 print("Downloading and extracting BLAST 2.7.1")
 os.system("wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz")
 os.system("tar -xzf ncbi-blast-2.7.1+-x64-linux.tar.gz && rm ncbi-blast-2.7.1+-x64-linux.tar.gz")
 blastn = os.getcwd() + "/ncbi-blast-2.7.1+/bin/blastn"

else:
 blastn = "blastn"

os.system("rm vers.txt")

if not os.path.exists("good-mapping"):
    os.makedirs("good-mapping")

if not os.path.exists("bad-mapping"):
    os.makedirs("bad-mapping")

if not os.path.exists("collapsed"):
    os.makedirs("collapsed")

#check if bowtie version 1.1.2 exist only of args.index is provided
if (args.mapping_stats == "YES") and (args.index != None):
   os.system("bowtie --version > vers.txt")
   infile= open("vers.txt", "r")
   lines = infile.readlines()
   if lines != []:
      ver1 = lines[0].strip().split("bowtie version")[1].split(" ")[1]
   else:
      ver1 = "no"
   if (ver1 != "1.1.2"):
      print("Downloading and extracting bowtie 1.1.2")
      os.system("wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip")
      os.system("unzip bowtie-1.1.2-linux-x86_64.zip && rm bowtie-1.1.2-linux-x86_64.zip")
      bowtie =  os.getcwd() + "/bowtie-1.1.2/bowtie"
      os.system("rm vers.txt")
   else:
      bowtie = "bowtie"



def worker(f):
            filename = f.split(".")[0]
            command = "mkdir "+ "-p "+ "aux_files/" + filename
            os.system(command)
            newpath = "aux_files/" + filename + "/"
            print("processing file " + filename + ".fastq") 
            log = newpath + filename + "_bowtie.txt"
            log2 = newpath + filename + "_cutadapt.txt"
            blast_file = newpath + filename + "_blast.csv"
            fastq_filename = filename + "_trimmed2.fastq"
            infile= open(f)
            fastq_lst = infile.readlines()[1::4]
            fastq_lst = [line.strip() for line in fastq_lst]
            collapsed = Counter(fastq_lst)
            collapsed = sorted(list(collapsed.items()), key=operator.itemgetter(1), reverse=True)
            if len(collapsed) == 0:
               sys.exit('\nERROR: Please check format of FASTQ files \n%s' %(docstring))
            fasta_file = "collapsed/" + filename + ".fa"
            ofile = open(fasta_file, "w")
            print("total number of unique sequences is " + str(len(collapsed)))
            ind_count = 0
            for x,y in collapsed:
              if (ind_count <= 700):
                 if (len(x) >= 29):
                  ind_count = ind_count +  1
                  ofile.write(">" +  str(y)+ ":Length=" + str(len(x))  + "|" + str(x) +  "\n" + str(x) + "\n")
                 else:
                  continue
              else:
                break
            ofile.close()
            
            command = "blastn -task blastn-short -ungapped -max_hsps 1 -max_target_seqs 1 -query " + fasta_file + " -subject " +  subject_file + " -outfmt '10 qseqid qseq length evalue qstart qend' -out " + blast_file + " &>dummy.out"
            print(command)
            os.system(command)
            csv = pd.read_csv(blast_file, header = None)
            csv.columns = ["sequence", "aligned_seq", "aligned_length", "evalue", "qstart", "qend"]
            csv = csv.loc[csv["aligned_length"] >= 21]
            if (len(csv) == 0):
               sys.exit('\nERROR: Check if the FASTQ files are already random-mer trimmed  \n%s' %(docstring))
            
            csv.insert(1,"seq_length", "")
            csv["sequence"] = csv["sequence"].apply(lambda x: x.split("|")[1])
            csv["seq_length"] = csv["sequence"].apply(lambda x: len(x))
            csv.insert(7,"rev_mer", "")
            csv["rev_mer"] = csv["seq_length"]-csv["qend"]        
            collapsed = Counter(csv["qstart"].tolist())
            collapsed = sorted(list(collapsed.items()), key=operator.itemgetter(1), reverse=True)
            collapsed = collapsed[:2]
            fivep_mer = int(collapsed[0][0]-1)
            if fivep_mer != 0:
             csv = csv.loc[csv["evalue"] >= 0.00014]
             #fivep_mer = csv["qstart"].mode()[0]
             print("processing file " + filename + ".fastq") 
             collapsed = Counter(csv["qstart"].tolist())
             collapsed = sorted(list(collapsed.items()), key=operator.itemgetter(1), reverse=True)
             collapsed = collapsed[:2]
             fivep_mer = int(collapsed[0][0]-1)
             collapsed = Counter(csv["rev_mer"].tolist())
             collapsed = sorted(list(collapsed.items()), key=operator.itemgetter(1), reverse=True)
             collapsed = collapsed[:2]
             threep_mer = int(collapsed[0][0])
            else:
             collapsed = Counter(csv["rev_mer"].tolist())
             collapsed = sorted(list(collapsed.items()), key=operator.itemgetter(1), reverse=True)
             collapsed = collapsed[:2]
             collapsed = sorted(collapsed, key=operator.itemgetter(0), reverse=True)
             threep_mer = int(collapsed[0][0])
            csv.to_csv(blast_file, sep = ",", index=False)
            if (fivep_mer != 0) and (threep_mer != 0):
             #if fivep_mer == threep_mer:
               
               command = "cutadapt -u " + str(fivep_mer) + " -u -" + str(threep_mer) + " -q 20 -m 15 -M50 " + " -o good-mapping/" + fastq_filename +  " " + f + " 2>> "+ log2
               print(command)
               os.system(command)
               if args.index != None:
                command = "bowtie --best -v 1 -p 20 " + index + " -q good-mapping/" + fastq_filename + " > check.sam 2> "+ log
                print(command)
                os.system(command)
                infile= open(log, "r")
                lines = infile.readlines()
                a = float(lines[1].strip().split("(")[1].split("%")[0])
                no_reads = lines[0].strip().split("processed: ")[1]
                if (a >=50):
                    element = [filename,fivep_mer,threep_mer,a,"good",no_reads]
                    return element
                else:
                    command = "mv good-mapping/"+ fastq_filename + " bad-mapping/" + fastq_filename
                    print(command)
                    os.system(command)
                    element = [filename,fivep_mer,threep_mer,a,"bad", no_reads]
                    return element
               else:
                element = [filename,fivep_mer,threep_mer,"na","na", "na"]
                return element
            elif (fivep_mer != 0) and (threep_mer == 0):
               command = "cutadapt -u " + str(fivep_mer) + " -q 20 -m 15 -M50 " + " -o good-mapping/" + fastq_filename +  " " + f + " 2>> "+ log2
               print(command)
               os.system(command)
               if args.index != None:
                command = "bowtie --best -v 1 -p 20 " + index + " -q good-mapping/" + fastq_filename + " > check.sam 2> "+ log
                print(command)
                os.system(command)
                infile= open(log, "r")
                lines = infile.readlines()
                a = float(lines[1].strip().split("(")[1].split("%")[0])
                no_reads = lines[0].strip().split("processed: ")[1]
                if (a >=50):
                    element = [filename,fivep_mer,"na",a,"good",no_reads]
                    return element
                else:
                    command = "mv good-mapping/"+ fastq_filename + " bad-mapping/" + fastq_filename
                    print(command)
                    os.system(command)
                    element = [filename,fivep_mer,"na",a,"bad", no_reads]
                    return element
               else:
                element = [filename,fivep_mer,"na","na","na", "na"]
                return element
            elif (fivep_mer == 0) and (threep_mer != 0):
               command = "cutadapt -u -" + str(threep_mer) + " -q 20 -m 15 -M50 " + " -o good-mapping/" + fastq_filename +  " " + f + " 2>> "+ log2
               print(command)
               os.system(command)
               if args.index != None:
                command = "bowtie --best -v 1 -p 20 " + index + " -q good-mapping/" + fastq_filename + " > check.sam 2> "+ log
                print(command)
                os.system(command)
                infile= open(log, "r")
                lines = infile.readlines()
                a = float(lines[1].strip().split("(")[1].split("%")[0])
                no_reads = lines[0].strip().split("processed: ")[1]
                if (a >=50):
                    element = [filename,"na",threep_mer,a,"good",no_reads]
                    return element
                else:
                    command = "mv good-mapping/"+ fastq_filename + " bad-mapping/" + fastq_filename
                    print(command)
                    os.system(command)
                    element = [filename,"na",threep_mer,a,"bad", no_reads]
                    return element
               else:
                element = [filename,"na",threep_mer,a,"na", "na"]
                return element

result_list =[]
pool = multiprocessing.Pool(processes = 20)
print("Processing " + str(len(files)) + " files")
result_list.append(pool.map(worker, [f for f in files if f.endswith(".fastq")]))
result_list = list(itertools.chain.from_iterable(result_list))
pool.close()
os.system("rm dummy.out")
asan  = pd.DataFrame.from_records(result_list, columns=['filename', '5p-randomer','3p-randomer','mapping-%','status','#reads'])       
asan.to_csv("random_mers.csv", sep = ",", index=False)
if args.index != None:
 command = "rm check.sam" 
 print(command)
 os.system(command)
            

