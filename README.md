# adapt_find
adapt_find identifies adapters sequences from single end raw sequencing files in FASTQ format. To run the script on raw FASTQ files, the following dependencies are required: PYTHON pandas module, scipy, numpy, cutadapt, blast, and bowtie.

# adapt_find usage 

adapt_find.py <sequencing platform> [-- min_len] [-- max_len] [-- index] [-- input_path] [-- output_path] [-- files]
  
"sequencing platform" is a mandatory argument and has to be specified. A list of the arguments and available options can be found [here](https://github.com/asan-nasa/adapt_find/blob/master/manual/adapt_find_manual.pdf) 

# adapt_find usage examples

For Illumina as sequencing platform and if the current directory has all FASTQ files

```$ python adapt_find.py ILLUMINA```

For Illumina as sequencing platform and to specify filenames explicitly

```$ python adapt_find.py ILLUMINA --files filename1.fastq filename2.fastq```

For Illumina as sequencing platform and to specify path to folder containing FASTQ files

```$ python adapt_find.py ILLUMINA --input_path path/to/folder```



# random_mer

random_mer identifies random_mer sequences from single-end adapter-trimmed FASTQ files. To run the script on raw FASTQ files, the following dependencies are required: PYTHON pandas module, cutadapt, blast, and bowtie.

# random_mer usage

```$ python random_mer.py path/to/genome/file.fa [--index path/to/index][--input_path path/to/folder][--output_path path/to/folder] [--files list of files]```

Example

For running all fastq files in current working directory

```$ python random_mer.py ~/genome/file.fa ```

For running only selected files

```$ python random_mer.py ~/genome/file.fa --files filename1.fastq filename2.fastq```


# fastqc_parser

fastqc_parser can parse the output of FASTQC reports: consolidate the results and extract the relevant information from all the output files. Depending on input arguments, fastqc_parser can either run fastqc on raw fastq files and then parse the output. Alternatively, it can be run directly on fastqc output files. To run this script, the following dependencies are required: PYTHON, pandas module, and FASTQC.

# fastqc_usage

```$ python fastqc_parser.py --run_fastqc [--fastqc_folder] [-- index]  [--param] [-- input_path] [-- output_path] [-- files] ```

Example

Run fastqc and parse the fastqc reports and store ouput in current working directory

```$ python fastqc_parser.py```


The different arguments and options for the three tools are described in the [manual](https://github.com/asan-nasa/adapt_find/tree/master/manual)




