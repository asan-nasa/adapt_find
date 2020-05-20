# adapt_find
adapt_find identifies adapters sequences from single-end raw sequencing files in FASTQ format. To run the script on raw FASTQ files, the following dependencies are required: PYTHON pandas module, scipy, numpy, cutadapt, blast, and bowtie.

# adapt_find usage 

adapt_find.py <sequencing_platform> [-- min_len] [-- max_len] [-- index] [-- input_path] [-- output_path] [-- files]
  
"sequencing_platform" is a mandatory argument and has to be specified. A list of arguments and available options can be found [here](https://github.com/asan-nasa/adapt_find/blob/master/manual/adapt_find_manual.pdf) 

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

For running random_mer on all FASTQ files in the current working directory

```$ python random_mer.py ~/genome/file.fa ```

For running random_mer only on selected FASTQ files

```$ python random_mer.py ~/genome/file.fa --files filename1.fastq filename2.fastq```


# fastqc_parser

fastqc_parser can parse the output of FASTQC reports: consolidate the results and extract the relevant information from all the output files. Depending on input arguments, fastqc_parser can either run fastqc on raw FASTQ files and then parse the output. Alternatively, it can be run directly on FASTQC output files. To run this script, the following dependencies are required: PYTHON, pandas module, and FASTQC.

# fastqc_parser usage

```$ python fastqc_parser.py --run_fastqc [--fastqc_folder] [-- index]  [--param] [-- input_path] [-- output_path] [-- files] ```

Example

Run FASTQC and parse the FASTQC reports. And store the output in the current working directory

```$ python fastqc_parser.py```


The different arguments and options for the three tools are described in the [manual](https://github.com/asan-nasa/adapt_find/tree/master/manual)




