# MITOPORE command line tool manual

### Getting start
Mitopore_workflow requires python3 and works on Unix liked environment
Simply install Python3 and run the code below(you may have to run it as a sudo user on Linux):
```console
git clone git@github.com:thachnguyen/mitopore_workflow.git
cd mitopore_workflow/mitopore_local
python install.py
```

#### Run mitopore_workflow pipeline
##### Data Preparation
Mitopore workflow supports multiple fastq files. The user must store all the fastq files in a single fastq folder.<br>
testdata_directory/fastq/sample1.fastq<br>
testdata_directory/fastq/sample2.fastq<br>
testdata_directory/fastq/sample3.fastq<br>

mitopore_workflow has two pipelines 
##### SNV calling
```console
python mitopore_snv.py testdata_directory_path 
```
##### INDEL calling (Beta version)
```console
python mitopore_indel.py testdata_directory_path 
```
##### Optional parameters
Running parameters are preset in config.yaml. 
##### Results
The result is summarized in a single HTML file report.html in your testdata_directory. Other supplementary result files (graphical plots, BAM alignment files, coverage, and mapping reports ...)are in the Results folder and Analysis folder.


### System requirement
    * CPU: 2.0 GHz (64bits) 2 cores or higher
    * Memory: 24 GB or higher
    * Diskdrive: 100 GB free space 
    * Linux (64 bits) or MacOS.
### Dependencies
All dependent packages are described in install.py files