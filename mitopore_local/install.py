import os
import shutil

#install samtools
os.system('sudo apt install samtools')
os.system('sudo apt install bcftools')
os.system('sudo apt install curl')
os.system('sudo apt install python3-pip')
os.system('sudo apt install default-jre')
os.system('pip install -r requirements.txt')

# download reference 
if not os.path.exists('./reference'):
    os.mkdir('reference')
print('Downloading reference genome')
os.system('wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz')
os.system('wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.MT.fa.gz')
os.system('wget https://ftp.ensembl.org/pub/release-110/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.MT.fa.gz')
os.system('wget https://ftp.ensembl.org/pub/release-110/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.chromosome.MT.fa.gz')
os.system('wget https://ftp.ensembl.org/pub/release-110/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.chromosome.MtDNA.fa.gz')
os.system('mv *.gz reference')
os.system('gunzip reference/*.gz')
print('All reference downloaded successfully')
if not os.path.exists('reference/bwa'):
    os.mkdir('reference/bwa')
    shutil.copyfile('reference/Homo_sapiens.GRCh38.dna.chromosome.MT.fa','reference/bwa/Homo_sapiens.GRCh38.dna.chromosome.MT.fa' )
if not os.path.exists('./tools'):
    os.mkdir('tools')
print('installing gatk')
os.chdir('tools')
os.system('wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip')
os.system('unzip gatk-4.3.0.0.zip')

print('Installing mutserve')
os.system('curl -sL mutserve.vercel.app | bash')

print('installing minimap2')
os.system('git clone https://github.com/lh3/minimap2')
os.chdir('minimap2')
os.system('make')
os.chdir('..')

print('installing bwa')
os.system('git clone https://github.com/lh3/bwa.git')
os.chdir('bwa')
os.system('make')
os.chdir('..')

os.system('wget https://github.com/genepi/haplogrep3/releases/download/v3.2.1/haplogrep3-3.2.1-linux.zip')
os.system('unzip haplogrep3-3.2.1-linux.zip')
os.chdir('..')

print('bwa indexing reference')
os.system('./tools/bwa/bwa index reference/bwa/Homo_sapiens.GRCh38.dna.chromosome.MT.fa')
os.system('samtools faidx reference/bwa/Homo_sapiens.GRCh38.dna.chromosome.MT.fa')
os.system('./tools/gatk-4.3.0.0/gatk CreateSequenceDictionary -R reference/bwa/Homo_sapiens.GRCh38.dna.chromosome.MT.fa')


# os.system('')
# os.system('pip install ')