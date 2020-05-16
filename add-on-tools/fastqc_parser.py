#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import pandas as pd
import multiprocessing
import argparse

docstring= """

USAGE:
python fastqc_parser.py  

Arguments:
Valid arguments:
(1) --fastqc_run = "YES" or "NO" (default: "YES")
(2) --fastqc_folder = PATH to fastqc output (applicable only if --fastqc_run = "YES")
(3) --input_path = PATH to FASTQ files (applicable only if --fastqc_run = "NO" and default current directory)
(4) --output_path = PATH to output folder to stored parsed output (If specified PATH doesn't exist, it will be created) default -folder named "output" in current directory)
(5) --param = either one or more of the following options: "PASS", "WARN", "FAIL" (Default: "FAIL")
(6) --files = enter FASTQ files seperated by space(if files are not in current directory,enter the absolute path and not the relative path)
Use --help for more info

Argument options for --fastqc_run and --param are case sensitive

DESCRIPTION
This tool parses and consolidates fastqc output

"""

parser = argparse.ArgumentParser(usage = "\n" +"python %(prog)s [--input_path path/to/folder] [--output_path path/to/folder] [--files list of files] [--param PASS,WARN,FAIL] \n" +"\n" + "\n\n" + "\n" + "Description:\nRun fastqc and parse its output\n", add_help=False)
optional = parser.add_argument_group('optional arguments')
optional.add_argument('--fastqc_run', default = "YES", help= 'enter YES if fastqc needs to be run, NO if the fastqc output is already there')
optional.add_argument('--fastqc_folder', help= 'paste path to input fastqc output files')
optional.add_argument('--input_path', default = os.getcwd(), help= 'paste path to input FASTQ files')
optional.add_argument('--output_path', default = os.getcwd()+"/output", help= 'paste path to store output files')
optional.add_argument('--param', nargs='*', default = ["FAIL"], help= 'enter one or more from the following [PASS,WARN,FAIL] ')
optional.add_argument('--files', nargs='*', help= 'enter FASTQ files seperated by space')
optional.add_argument("-h", "--help", action='help', help='print help message')
args = parser.parse_args()

if (args.fastqc_run != "YES") and (args.fastqc_folder != None):
  folder = args.fastqc_folder
else:
  folder = args.output_path
  if not os.path.exists(folder):
    os.makedirs(folder)
  if not os.path.exists(folder+"/fastqc_reports"):
    os.makedirs(folder+"/fastqc_reports")


filter_lst = args.param

if (args.fastqc_run != "YES"): 
 if (args.fastqc_run != "NO"):
   sys.exit('\nERROR: Argument for fastqc_run is not valid: Allowed arguments for fastqc_run is "YES" or "NO" \n%s' %(docstring))

if (args.fastqc_run != "YES") and (args.fastqc_folder == None):
   sys.exit('\nERROR: Argument for fastqc_run is not valid. --fastqc_folder option with path to fastqc reports have to specified \n%s' %(docstring))

if (args.fastqc_run == "YES") and (args.fastqc_folder != None):
   sys.exit('\nERROR: Argument for fastqc_run is not valid. if --fastqc_folder option is specified. argument fastqc_run should be specified as "YES"  \n%s' %(docstring))
   

for elem in filter_lst:
    if elem in ["FAIL", "WARN", "PASS"]:
       continue
    else:
       print("Invalid Argument: Check case, spelling or available arugument options" + elem)
       sys.exit('\nERROR: Argument for fastqc_run is not valid\n%s' %(docstring))
        

if (args.files == None) and (args.fastqc_run == "YES"):
 cwd = args.input_path
 files = [f for f in os.listdir(cwd) if f.endswith(".fastq")]
 file_no = len(files)
elif (args.files != None) and (args.fastqc_run == "YES"): 
 files = args.files
 file_no = len(files)



if (args.fastqc_run == "YES") and (len(files) == 0):
  sys.exit('\nERROR: No fastq files found  \n%s' %(docstring))

def worker(f):
  filename = f.split(".")[0]
  print("processing file " + filename + ".fastq") 
  command = "fastqc " + f + " -o " + folder+"/fastqc_reports"
  os.system(command)

if args.fastqc_run == "YES":
   pool = multiprocessing.Pool(processes = file_no)
   pool.map(worker, [f for f in files]) 
   pool.close()
   pool.join()
elif args.fastqc_run == "NO":
   if args.fastqc_folder == None:
      sys.exit('\nERROR: Argument for fastqc_run is not valid\n%s' %(docstring))
else:
   sys.exit('\nERROR: Argument for fastqc_run is not valid\n%s' %(docstring))

os.chdir(folder+"/fastqc_reports")
files2 = [f for f in os.listdir(os.getcwd()) if f.endswith(".zip")]
if len(files2) == 0:
  sys.exit('\nERROR: No fastqc output zip files found. Please check if the specified path to fastqc reports folder is correct' %(docstring)) 
os.system("unzip \*.zip")
os.chdir("../..")


if not os.path.exists(args.output_path):
    os.makedirs(args.output_path)

print(folder+"/fastqc_reports")
master = []
master_dict = {"Adapter Content":"adapter_content.png","Sequence Duplication Levels":"duplication_levels.png","Per base N content":"per_base_n_content.png","Per base sequence content":"per_base_sequence_content.png","Per sequence quality scores":"per_sequence_quality.png","Per tile sequence quality":"per_tile_quality.png","Sequence Length Distribution":"sequence_length_distribution.png"}

for root, dirs, files in os.walk(folder):
     for name in files:
         if name == "summary.txt":
            ind = []
            filepath = os.path.join(root, name)
            df = pd.read_csv(filepath,sep = '\t', header = None)
            df.columns = ["result", "parameters", "filename"]
            #print df["parameters"].tolist()
            ind.append(df["filename"].tolist()[0])
            ind.extend(df["result"].tolist())
            master.append(ind)
            #if len(sys.arg.param) 
            for par in filter_lst:
             problem = df.loc[df["result"]==par]
             awk_path = filepath.split("summary.txt")[0] + "fastqc_data.txt" 
             of_name = filepath.split("summary.txt")[0].split("/")[-2]
             copypath =  filepath.split("summary.txt")[0] + "Images/"
             newfilepath = args.output_path + "/"+ par + "/"+ filepath.split("summary.txt")[0].split("/")[-2] + "/"
             module = ["Basic Statistics", "Per base sequence quality", "Per tile sequence quality", "Per sequence quality scores", "Per base sequence content", "Per sequence GC content", "Per base N content","Sequence Length Distribution","Sequence Duplication Levels", "Overrepresented sequences", "Adapter Content"]
             if len(problem)!=0:
              if not os.path.exists(args.output_path+ "/"+ par + "/"+ of_name):
                os.makedirs(args.output_path+ "/"+ par + "/"+ of_name)
              for item in (problem["parameters"].tolist()):
                if item in master_dict:
                   command = "cp " + copypath + master_dict[item] + " " + newfilepath + master_dict[item]
                   os.system(command)
                   command = "awk '/" + item +"/{flag=1; next} /END_MODULE/{flag=0} flag' " + awk_path + " > " + newfilepath + master_dict[item].split(".")[0] +  ".tsv"
                   #print command
                   os.system(command)
                else:
                   continue

if len(master) != 0:
 col_names = ["filename"]
 col_names.extend(df["parameters"].tolist())          
 asan  = pd.DataFrame.from_records(master, columns=col_names)
 result = asan.drop("filename", axis=1)
 result = result.apply(pd.value_counts).fillna(0)
 result.to_csv(args.output_path+"/ind_count.csv", sep = ",")
 asan.to_csv(args.output_path+ "/consolidated_results_fastqc.csv", sep = ",", index=False)

else:
 sys.exit('\nERROR: No extracted fastqc output files found. Check fastqc ouput files\n%s' %(docstring))
