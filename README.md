# MITOPORE commanline tool manual

### Getting start
Mitopore_workflow require python3 and work on Unix liked environment
Simply install Python3 and run code below(you may have to run it as sudo user on Linux):
```console
gitclone git@github.com:thachnguyen/mitopore_workflow.git
cd mitopore_workflow/mitopore_local
python install.py
```

#### Run mitopore_workflow pipeline
##### Data preparation
Mitopore workflow support multiple fastq file. User must store all the fastq file into a single fastq folder.<br>
testdata_directory/fastq/sample1.fastq<br>
testdata_directory/fastq/sample2.fastq<br>
testdata_directory/fastq/sample3.fastq<br>

mitopore_workflow has two pipelines 
##### SNV calling
```console
python mitopore_snv.py testdata_directory_path 
```
##### INDEL calling
```console
python mitopore_indel.py testdata_directory_path 
```
### System requirement
#### System requirement
    * CPU: 2.0 GHz (64bits) 2 cores or higher
    * Memory: 24 GB or higher
    * Diskdrive: 100 GB free space 
    * Linux (64 bits) or MacOS.

### Dependancies
All dependancies package are described in install.py files