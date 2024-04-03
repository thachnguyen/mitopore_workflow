'''
This is the main application functions. Read the user's upload files (fastaq format).'
'''
import os
import pandas as pd
import json
import shutil
import yaml
from Bio import SeqIO

def create_yaml(path, samples, yaml_file = 'config.yaml'):
    # '''
    # Create the user defined YAML from YAML template for RScript. The fastq files are separate into different groups as subfolders. Named by directory's name and file'sname.
    # Create auxilary file for mutserver
    # '''   
    with open(yaml_file) as file:
        default_config = yaml.load(file, Loader=yaml.FullLoader)
    # with open('static/config.yaml') as file1:
    #     default_config1 = yaml.load(file1, Loader=yaml.FullLoader)

    default_config['Samples'] = samples
    with open('%s/config.yaml'%path, 'w') as f:
        yaml.dump(default_config, f)
    return

def manage_fastq_list(path):
    '''
    Read user upload files in zip format and return the list of fastq entries for later process, pass it in to YAML file
    Run Quality control
    If BAM files uploaded, skip alignment step
    '''
    samples_data = os.listdir(path+'/fastq')
    # accept fastq and fastq.gz, bams only
    sample_1 = []
    for sample in samples_data:
        if sample[-5:] == 'fastq':
            sample_1.append(sample)
        elif sample[-8:]== 'fastq.gz':
            os.system('gunzip %s/fastq/%s'%(path,sample))
            sample_1.append(sample[:-3])
        elif sample[-3:]== 'bam':
            sample_1.append(sample)
        else:
            if os.path.isfile('%s/fastq/%s'%(path, sample)):
                os.remove('%s/fastq/%s'%(path, sample))
            else:
                shutil.rmtree('%s/fastq/%s'%(path, sample))
    samples_data = sample_1
    group1 = {}
    for i, each_file in enumerate(samples_data):
        group1[each_file.split('.')[0]] = 'fastq/%s'%each_file

    return group1

def run_minimap2(path='data', organism = 'human', seq_data = 'nanopore'):
    file_org={'human':'Homo_sapiens.GRCh38.dna.chromosome.MT.fa',
                'rat': 'Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.MT.fa',
                'mouse':'Mus_musculus.GRCm39.dna.chromosome.MT.fa',
                'zebrafish':'Danio_rerio.GRCz11.dna.chromosome.MT.fa',
                'celegans':'Caenorhabditis_elegans.WBcel235.dna.chromosome.MtDNA.fa',
                'custom': 'custom_ref.fasta'}
    if not os.path.exists('%s/fastq'%path):
        os.mkdir('%s/fastq'%path)
    os.system('mv %s/*.fastq %s/fastq'%(path, path))
    if os.listdir('%s/fastq'%path) == 0:
        print('No fastq files, please check your fastq data')
    if not os.path.exists('%s/Analysis'%path):
        os.mkdir('%s/Analysis'%path)
    path_minimap = '%s/Analysis/Minimap'%path
    if not os.path.exists(path_minimap):
        os.mkdir(path_minimap)

    path_minimap_shifted = '%s/Analysis/Minimap/shifted'%path
    if not os.path.exists(path_minimap_shifted):
        os.mkdir(path_minimap_shifted)

    path_flagstat = '%s/Analysis/flagstat/'%path
    if not os.path.exists(path_flagstat):
        os.mkdir(path_flagstat)
    
    path_depths = '%s/Analysis/depths/'%path
    if not os.path.exists(path_depths):
        os.mkdir(path_depths)
    with open('%s/Analysis/sample.txt'%path, 'w') as f1:
        f1.write('')
    for fastq_file in os.listdir('%s/fastq/'%(path)):
        path2 = '%s/fastq/%s'%(path, fastq_file)
        fastq_file1 = fastq_file.split('.')[0]
        with open('%s/Analysis/sample.txt'%path, 'a') as f1:
            f1.write(fastq_file1+'.bam\n')
        if seq_data == 'nanopore':
            os.system('tools/minimap2/minimap2 -t 4 -ax map-ont --secondary=no ./reference/%s %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
            print('minimap2 -t 4 -ax map-ont --secondary=no *** ./reference/%s %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
            if organism == 'human':
                os.system('minimap2 -t 4 -ax map-ont --secondary=no /home/ag-rossi/projects/mitopore_workflow/mitopore_local/reference/mtDNA/rCRS_shifted.fasta %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(path2, path_minimap_shifted, fastq_file1))
        
        elif seq_data == 'pacbio':
            os.system('tools/minimap2/minimap2 -t 4 -ax map-pb --secondary=no ./reference/%s %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
            print('minimap2 -t 4 -ax map-pb -uf --secondary=no *** ./reference/%s %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
            if organism == 'human':
                os.system('minimap2 -t 4 -ax map-pb -uf --secondary=no /home/ag-rossi/projects/mitopore_workflow/mitopore_local/reference/mtDNA/rCRS_shifted.fasta %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(path2, path_minimap_shifted, fastq_file1))

        elif seq_data =='illumina':
            os.system('bwa mem -K 100000000 -p -v 3 -t 4 -Y ./reference/%s %s |samtools view -Sb|samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
            print('bwa mem -K 100000000 -p -v 3 -t 4 -Y ./reference/%s %s |samtools view -Sb|samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
            if organism == 'human':
                os.system('bwa mem -K 100000000 -p -v 3 -t 4 -Y /home/ag-rossi/projects/mitopore_workflow/mitopore_local/reference/mtDNA/rCRS_shifted.fasta %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(path2, path_minimap_shifted, fastq_file1))        
        else:
            print('sequence method is not listed')
        os.system('samtools flagstat %s/%s.bam>%s%s.txt'%(path_minimap, fastq_file1, path_flagstat, fastq_file1))
        # Calculating depth tracks
        print('Run depth calculation for %s'%fastq_file) 
        os.system('samtools depth %s/%s.bam>%s%s.txt'%(path_minimap, fastq_file1, path_depths, fastq_file1))
        # dump depth info to json, for CGViewJS part 
    return

def run_mutect2(path='data', organism = 'human'):
    file_org={'human':'Homo_sapiens.GRCh38.dna.chromosome.MT.fa',
                'rat': 'Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.MT.fa',
                'mouse':'Mus_musculus.GRCm39.dna.chromosome.MT.fa',
                'zebrafish':'Danio_rerio.GRCz11.dna.chromosome.MT.fa',
                'celegans':'Caenorhabditis_elegans.WBcel235.dna.chromosome.MtDNA.fa',
                'custom': 'custom_ref.fasta'}

    if os.listdir('%s/fastq'%path) == 0:
        print('No fastq files, please check your fastq data')
    if not os.path.exists('%s/Analysis'%path):
        os.mkdir('%s/Analysis'%path)
    path_minimap = '%s/Analysis/Minimap'%path
    if not os.path.exists(path_minimap):
        os.mkdir(path_minimap)
    path_minimap_shifted = '%s/Analysis/Minimap/shifted'%path
    if not os.path.exists(path_minimap_shifted):
        os.mkdir(path_minimap_shifted)

    path_flagstat = '%s/Analysis/flagstat/'%path
    if not os.path.exists(path_flagstat):
        os.mkdir(path_flagstat)
    
    path_depths = '%s/Analysis/depths/'%path
    if not os.path.exists(path_depths):
        os.mkdir(path_depths)

    path_mutect = '%s/Analysis/INDEL/'%path
    if not os.path.exists(path_mutect):
        os.mkdir(path_mutect)

    for fastq_file in os.listdir('%s/fastq_fil/'%(path)):
        path2 = '%s/fastq_fil/%s'%(path, fastq_file)
        fastq_file1 = fastq_file.split('.')[0]
        with open('%s/Analysis/sample.txt'%path, 'a') as f1:
            f1.write(fastq_file1+'.bam\n')
        # This protocol, we all use bwa-mem as aligner 
        os.system('./tools/bwa/bwa mem -K 100000000 -p -v 3 -t 4 -Y reference/bwa/%s %s|samtools view -Sb|samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
        os.system('./tools/gatk-4.3.0.0/gatk AddOrReplaceReadGroups -I %s/%s.bam -O %s/%s_addedRG.bam --RGID GT19-38445 --RGLB wgs --RGPL illumina --RGPU GT19-38445 --SORT_ORDER coordinate --CREATE_INDEX true --TMP_DIR /tmp/ --RGSM wgs1'%(path_minimap,fastq_file1, path_minimap, fastq_file1))
        os.system('./tools/gatk-4.3.0.0/gatk Mutect2 -R reference/bwa/Homo_sapiens.GRCh38.dna.chromosome.MT.fa --mitochondria-mode -I %s/%s_addedRG.bam -O %s/INDEL_%s.vcf'%(path_minimap, fastq_file1, path_minimap,fastq_file1))
        os.unlink('%s/%s_addedRG.bam'%(path_minimap, fastq_file1))
        os.unlink('%s/%s_addedRG.bai'%(path_minimap, fastq_file1))
        shutil.copyfile('%s/INDEL_%s.vcf'%(path_minimap,fastq_file1), '%sINDEL_%s.vcf'%(path_mutect,fastq_file1))
        os.system('samtools flagstat %s/%s.bam>%s%s.txt'%(path_minimap, fastq_file1, path_flagstat, fastq_file1))
    # Calculating depth tracks
        print('Run depth calculation for %s'%fastq_file) 
        os.system('samtools depth %s/%s.bam>%s%s.txt'%(path_minimap, fastq_file1, path_depths, fastq_file1))
        # dump depth info to json, for CGViewJS part 
    return

def write_rscript(path='data'):
    current_dir = os.getcwd()
    new_R = 'setwd("%s")\nref_dir<-"%s"\n'%(path, current_dir)
    new_R += open('R/RNA1.R', 'r').read()
    f = open(path+'/Analysis/RNA.R', 'w')
    f.write(new_R)
    f.close()
    return
        
def run_mutserver(path = 'data',organism = 'human', thres = '0.05'):
    file_org={'human':'rCRS.fasta',
            'rat': 'Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.MT.fa',
            'mouse':'Mus_musculus.GRCm39.dna.chromosome.MT.fa',
            'zebrafish':'Danio_rerio.GRCz11.dna.chromosome.MT.fa',
            'celegans':'Caenorhabditis_elegans.WBcel235.dna.chromosome.MtDNA.fa',
            'custom': 'custom_ref.fasta'}
    if not os.path.exists('%s/Analysis/Results/'%path):
        os.mkdir('%s/Analysis/Results/'%path)
    os.system('./tools/mutserve call %s/Analysis/Minimap/*.bam --reference reference/%s --output result1.vcf --threads 4 --baseQ 10 --level %s'%(path, file_org[organism], thres))    
    
    df1 = pd.read_csv('result1.txt', delimiter='\t')
    df11 = pd.read_csv('result1.vcf', delimiter='\t', skiprows= 7)
    if organism=='human':
        os.system('./tools/mutserve call %s/Analysis/Minimap/shifted/*.bam --reference /home/ag-rossi/projects/mitopore_workflow/mitopore_local/reference/mtDNA/rCRS_shifted.fasta --output result1_shifted.vcf --threads 4 --baseQ 10 --level %s'%(path, thres))    
        df2 = pd.read_csv('result1_shifted.txt', delimiter='\t')
        df2['Pos'] = df2['Pos'].apply(lambda x: x - 8569 if x > 8569 else x + 8000)
        df1 = pd.concat([df1, df2])
        df1 = df1.drop_duplicates()
        df1.to_csv('result1.txt', sep='\t', index=False)

        df22 = pd.read_csv('result1_shifted.vcf', delimiter='\t', skiprows= 7)
        df22['POS'] = df22['POS'].apply(lambda x: x - 8569 if x > 8569 else x + 8000)
        df11 = pd.concat([df11, df22])
        df11 = df11.drop_duplicates()
        df11.to_csv('result1.vcf', sep='\t', index=False)
        header = '##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="Variants passed mtDNA-Server">\n##FORMAT=<ID=AF,Number=1,Type=String,Description="Inferred Allele Frequency of top (non-reference) allele">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##Mutserve=v2.0.0-rc13\n##contig=<ID=chrM,length=16569>\n'
        with open('result1.vcf', 'r') as file:
            existing_content = file.read()
        new_content = header + existing_content
        with open('result1.vcf', 'w') as file:
            file.write(new_content)

    os.system('mv result1.vcf %s/Analysis/Results/'%path)
    os.system('mv result1.txt %s/Analysis/Results/'%path)
    os.system('mv result1_shifted.vcf %s/Analysis/Results/'%path)
    os.system('mv result1_shifted.txt %s/Analysis/Results/'%path)

    os.system('bcftools view -S %s/Analysis/sample.txt %s/Analysis/Results/result1.vcf > %s/Analysis/Results/final_output_2percent.vcf --force-samples'%(path,path,path))
    df = pd.read_csv('%s/Analysis/Results/final_output_2percent.vcf'%path, delimiter='\t', skiprows=11)
    df = df.drop(['ID', 'QUAL', 'INFO'], axis=1)
    df.to_html('%s/Analysis/Results/VCF.html'%path)
    return

def read_indel(path='data'):
    indel_all = pd.DataFrame()
    for sample1 in os.listdir('%s/Analysis/INDEL'%(path)):
        indel_df = pd.read_csv('%s/Analysis/INDEL/%s'%(path, sample1), delimiter='\t')
        indel_df['ID'] = sample1.split('.')[0][6:]
        indel_all = pd.concat([indel_all, indel_df], ignore_index=True, sort=False)
    indel_all.to_html('%s/Analysis/Results/INDEL.html'%path)

    return indel_all

def write_json_cgview(path='data', baseline_cgv = 1000):
    # For json plots
    my_json = json.load(open('static/JS_library/NZ_mito.json'))
    plot_list = []
    track_list = my_json['cgview']['tracks']
    features_list = my_json['cgview']['features']
    path_depths = '%s/Analysis/depths/'%path
    vcf_list = pd.read_csv('%s/Analysis/Results/result1.txt'%path, delimiter='\t')

    for fastq_file in os.listdir('%s/fastq/'%(path)):
        fastq_file1 = fastq_file.split('.')[0]
        seq_depth = pd.read_csv('%s%s.txt'%(path_depths, fastq_file1), delimiter='\t', header=None)
        seq_depth_list = {}
        seq_depth_list['type'] = 'line'
        # Baseline of depth plots track
        seq_depth_list['baseline'] = baseline_cgv
        #Avoid filename bugs
        if len(str(fastq_file))>6:
            seq_depth_list['source'] = str(fastq_file)[:-6]
        else:
            seq_depth_list['source'] = str(fastq_file)

        seq_depth_list['legendPositive'] = 'High coverage'
        seq_depth_list['legendNegative'] = 'Low coverage'
        seq_depth_list['axisMin'] = 0
        seq_depth_list['positions'] = list(seq_depth[1])
        seq_depth_list['scores'] = list(seq_depth[2])
       
        plot_list.append(seq_depth_list)
        track1 = {'name': 'Coverage %s'%seq_depth_list['source'],
            'separateFeaturesBy': 'strand',
            'position': 'inside',
            'thicknessRatio': 1,
            'dataType': 'plot',
            'dataMethod': 'source',
            'dataKeys': '%s'%seq_depth_list['source']}
    
        track_list.append(track1)
        variance_list = vcf_list[(vcf_list['ID']=='%s.bam'%str(fastq_file1))&(vcf_list['Filter']=='PASS')]
        for i1 in range(len(variance_list)):
            feature1 = {'type': 'SNV',
                        'name': '%s %s>%s %.2f'%(fastq_file1, list(variance_list['Ref'])[i1], list(variance_list['Variant'])[i1],list(variance_list['VariantLevel'])[i1]),
                        'start': list(variance_list['Pos'])[i1],
                        'stop': list(variance_list['Pos'])[i1],
                        'strand': -1,
                        'source': str(fastq_file),
                        'contig': 'NC_012920',
                        'legend': 'SNV'}
            features_list.append(feature1)

        track2 = {'name': 'Variance %s'%seq_depth_list['source'],
            'separateFeaturesBy': 'strand',
            'position': 'inside',
            'thicknessRatio': 1,
            'dataType': 'plot',
            'dataMethod': 'source',
            'dataKeys': str(fastq_file)}
        track_list.append(track2)

    my_json['cgview']['plots'] = plot_list
    my_json['cgview']['tracks'] = track_list
    my_json['cgview']['features'] = features_list
    json_str = json.dumps(my_json)
    js_str = 'var NZ_mito = %s'%json_str
    with open('%s/Analysis/my_cgview.js'%path, 'w') as outfile:
        outfile.write(js_str)
    return

def write_json_cgview_indel(path='data', baseline_cgv = 1000):
    # For json plots
    my_json = json.load(open('static/JS_library/NZ_mito.json'))
    plot_list = []
    track_list = my_json['cgview']['tracks']
    features_list = my_json['cgview']['features']
    path_depths = '%s/Analysis/depths/'%path

    for fastq_file in os.listdir('%s/fastq/'%(path)):
        fastq_file1 = fastq_file.split('.')[0]
        seq_depth = pd.read_csv('%s%s.txt'%(path_depths, fastq_file1), delimiter='\t', header=None)
        seq_depth_list = {}
        seq_depth_list['type'] = 'line'
        # Baseline of depth plots track
        seq_depth_list['baseline'] = baseline_cgv
        #Avoid filename bugs
        if len(str(fastq_file))>6:
            seq_depth_list['source'] = str(fastq_file)[:-6]
        else:
            seq_depth_list['source'] = str(fastq_file)

        seq_depth_list['legendPositive'] = 'High coverage'
        seq_depth_list['legendNegative'] = 'Low coverage'
        seq_depth_list['axisMin'] = 0
        seq_depth_list['positions'] = list(seq_depth[1])
        seq_depth_list['scores'] = list(seq_depth[2])
       
        plot_list.append(seq_depth_list)
        track1 = {'name': 'Coverage %s'%seq_depth_list['source'],
            'separateFeaturesBy': 'strand',
            'position': 'inside',
            'thicknessRatio': 1,
            'dataType': 'plot',
            'dataMethod': 'source',
            'dataKeys': '%s'%seq_depth_list['source']}
    
        track_list.append(track1)
        indel_list = read_indel(path=path)
        for i1 in range(len(indel_list)):
            feature2 = {'type': 'INDEL',
                        'name': '%s %s>%s %.2f'%(list(indel_list['ID'])[i1], list(indel_list['REF'])[i1], list(indel_list['ALT'])[i1],list(indel_list['AF'])[i1]),
                        'start': list(indel_list['POS'])[i1],
                        'stop': list(indel_list['POS'])[i1],
                        'strand': -1,
                        'source': list(indel_list['ID'])[i1],
                        'contig': 'NC_012920',
                        'legend': 'INDEL'}
            features_list.append(feature2)

        track2 = {'name': 'Variance %s'%seq_depth_list['source'],
            'separateFeaturesBy': 'strand',
            'position': 'inside',
            'thicknessRatio': 1,
            'dataType': 'plot',
            'dataMethod': 'source',
            'dataKeys': str(fastq_file)}
        track_list.append(track2)

    my_json['cgview']['plots'] = plot_list
    my_json['cgview']['tracks'] = track_list
    my_json['cgview']['features'] = features_list
    json_str = json.dumps(my_json)
    js_str = 'var NZ_mito1 = %s'%json_str
    with open('%s/Analysis/my_cgview1.js'%path, 'w') as outfile:
        outfile.write(js_str)

    return


def select_subset(seq_record, min_length = 50, min_quality = 6):
    """_summary_ transform seq record into smaller segment which has higher quality than min_quality

    Args:
        seq_record (_type_): _description_ Seq record from Biopython
        min_length (int, optional): _description_. Defaults to 50.
        min_quality (int, optional): _description_. Defaults to 6.

    Returns:
        _type_: list 
        _description_ list of index 
    """
    store_list=[]
    phred_list = seq_record.letter_annotations['phred_quality']
    store_list_temp = []
    for i, base in enumerate(phred_list):
        if base > min_quality:
            store_list_temp.append(i)
        else:
            if len(store_list_temp)> min_length:
                store_list.append(store_list_temp)
                store_list_temp=[]
            else:
                store_list_temp=[]
    return store_list

def convert_fastq(fastqfile, outfile ='output.fastq'):
    """_summary_

    Args:
        fastqfile (_type_): input fastq to convert

    Returns:
        _type_: _description_
    """
    my_fastq = list(SeqIO.parse(fastqfile, format='fastq'))
    convert_record = []
    for seq_record in my_fastq:
        store_list = select_subset(seq_record)
        for i, item1 in enumerate(store_list):
            temp_seq = seq_record[item1[0]:item1[-1]]
            temp_seq.id = temp_seq.id + '_%s'%str(i)
            convert_record.append(temp_seq)
            # convert_record.append(seq_record[item1[0]:item1[-1]])
    for read in convert_record:
        with open(outfile, "a") as out_handle:
            SeqIO.write(read, out_handle, "fastq")

    return convert_record

def write_fastq(record_list, outfile = 'output.fastq'):
    """_summary_
    Args:
        record_list (_type_): _description_
        outfile (str, optional): _description_. Defaults to 'output.fastq'.
    """
    for read in record_list:
        with open(outfile, "a") as out_handle:
            SeqIO.write(read, out_handle, "fastq")
    return

def run_haplogrep3(path, tree = 'phylotree-rcrs@17.2', outfile = "haplo.txt"):
    os.system("./tools/haplogrep3 classify --in %s/Analysis/Results/result1.vcf --tree %s --out %s --extend-report"%(path, tree, outfile))
    df = pd.read_csv(outfile, delimiter='\t')
    df.to_html('%s/Analysis/Results/haplogrep.html'%path)
    return

