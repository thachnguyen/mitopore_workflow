import os

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


os.system('')
os.system('pip install ')