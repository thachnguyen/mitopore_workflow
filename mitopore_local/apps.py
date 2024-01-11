'''
This is the main application functions. Read the user's upload files (fastaq format).'
'''
import os
import pandas as pd
import json
import shutil

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
    if os.listdir('%s/fastq') == 0:
        print('No fastq files, please check your fastq data')
    if not os.path.exists('%s/Analysis'%path):
        os.mkdir('%s/Analysis'%path)
    path_minimap = '%s/Analysis/Minimap'%path
    if not os.path.exists(path_minimap):
        os.mkdir(path_minimap)

    path_flagstat = '%s/Analysis/flagstat/'%path
    if not os.path.exists(path_flagstat):
        os.mkdir(path_flagstat)
    
    path_depths = '%s/Analysis/depths/'%path
    if not os.path.exists(path_depths):
        os.mkdir(path_depths)

    for fastq_file in os.listdir('%s/fastq/'%(path)):
        path2 = '%s/fastq/%s'%(path, fastq_file)
        fastq_file1 = fastq_file.split('.')[0]
        if seq_data == 'nanopore':
            os.system('minimap2 -t 4 -ax map-ont --secondary=no ./reference/%s %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
            print('minimap2 -t 4 -ax map-ont --secondary=no *** ReferenceData/%s %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
        elif seq_data == 'pacbio':
            os.system('minimap2 -t 4 -ax map-pb --secondary=no ./reference/%s %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
            print('minimap2 -t 4 -ax map-pb -uf --secondary=no *** ReferenceData/%s %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
        elif seq_data =='illumina':
            os.system('bwa mem -K 100000000 -p -v 3 -t 4 -Y ./reference/%s %s |samtools view -Sb|samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
            print('bwa mem -K 100000000 -p -v 3 -t 4 -Y ./reference/%s %s |samtools view -Sb|samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
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

    if not os.path.exists('%s/fastq'%path):
        os.mkdir('%s/fastq'%path)
    os.system('mv %s/*.fastq %s/fastq'%(path, path))
    if os.listdir('%s/fastq') == 0:
        print('No fastq files, please check your fastq data')
    if not os.path.exists('%s/Analysis'%path):
        os.mkdir('%s/Analysis'%path)
    path_minimap = '%s/Analysis/Minimap'%path
    if not os.path.exists(path_minimap):
        os.mkdir(path_minimap)

    path_flagstat = '%s/Analysis/flagstat/'%path
    if not os.path.exists(path_flagstat):
        os.mkdir(path_flagstat)
    
    path_depths = '%s/Analysis/depths/'%path
    if not os.path.exists(path_depths):
        os.mkdir(path_depths)

    path_mutect = '%s/Analysis/INDEL/'%path
    if not os.path.exists(path_mutect):
        os.mkdir(path_mutect)

    for fastq_file in os.listdir('%s/fastq/'%(path)):
        path2 = '%s/fastq/%s'%(path, fastq_file)
        fastq_file1 = fastq_file.split('.')[0]
        # This protocol, we all use bwa-mem as aligner 
        os.system('bwa mem -K 100000000 -p -v 3 -t 4 -Y ./reference/bwa/%s %s|samtools view -Sb|samtools sort - -o %s/%s.bam'%(file_org[organism], path2, path_minimap, fastq_file1))
        os.system('./tools/GATK/gatk-4.3.0.0/gatk AddOrReplaceReadGroups -I %s/%s.bam -O %s/%s_addedRG.bam --RGID GT19-38445 --RGLB wgs --RGPL illumina --RGPU GT19-38445 --SORT_ORDER coordinate --CREATE_INDEX true --TMP_DIR /tmp/ --RGSM wgs1'%(path_minimap,fastq_file1, path_minimap, fastq_file1))
        os.system('./tools/GATK/gatk-4.3.0.0/gatk Mutect2 -R /home/ag-rossi/ReferenceData/Homo_sapiens.GRCh38.dna.chromosome.chrM.fa --mitochondria-mode -I %s/%s_addedRG.bam -O %s/INDEL_%s.vcf'%(path_minimap, fastq_file1, path_minimap,fastq_file1))
        os.unlink('%s/%s_addedRG.bam'%(path_minimap, fastq_file1))
        os.unlink('%s/%s_addedRG.bai'%(path_minimap, fastq_file1))
        shutil.copyfile('%s/INDEL_%s.vcf'%(path_minimap,fastq_file1), '%sINDEL_%s.vcf'%(path_mutect,fastq_file1))
        os.system('samtools flagstat %s/%s.bam>%s%s.txt'%(path_minimap, fastq_file1, path_flagstat, fastq_file1))
    # Calculating depth tracks
        print('Run depth calculation for %s'%fastq_file) 
        os.system('samtools depth %s/%s.bam>%s%s.txt'%(path_minimap, fastq_file1, path_depths, fastq_file1))
        # dump depth info to json, for CGViewJS part 
    return

def write_rscript(path='users_file/', s_id = 'Test_name_1618217069/'):
    new_R = 'setwd("/home/ag-rossi/projects/mitopore/mitopore/%s%s")\n'%(path, s_id)
    new_R += open('R/RNA1.R', 'r').read()
    f = open(path+s_id+'/Analysis/RNA.R', 'w')
    f.write(new_R)
    f.close()
    return
        
def run_mutserver(organism = 'human'):
    file_org={'human':'rCRS.fasta',
            'rat': 'Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.MT.fa',
            'mouse':'Mus_musculus.GRCm39.dna.chromosome.MT.fa',
            'zebrafish':'Danio_rerio.GRCz11.dna.chromosome.MT.fa',
            'celegans':'Caenorhabditis_elegans.WBcel235.dna.chromosome.MtDNA.fa',
            'custom': 'custom_ref_%s.fasta'%s_id}
    os.system('./mutserver/mutserve call users_file/%s/Analysis/Minimap/*.bam --reference /home/ag-rossi/projects/mitopore/mitopore/reference/%s --output users_file/%s/Analysis/Results/result1.vcf --threads 4 --baseQ 10 --level 0.05'%(s_id, file_org[organism] ,s_id))    
    os.system('bcftools view -S users_file/%s/Analysis/sample.txt users_file/%s/Analysis/Results/result1.vcf > users_file/%s/Analysis/Results/final_output_2percent.vcf --force-samples'%(s_id,s_id,s_id))
    df = pd.read_csv('users_file/%s/Analysis/Results/final_output_2percent.vcf'%s_id, delimiter='\t', skiprows=11)
    df = df.drop(['ID', 'QUAL', 'INFO'], axis=1)
    df.to_html('users_file/%s/Analysis/Results/VCF.html'%s_id)
    return

def read_indel(path='users_file/', s_id = 'Test_name_1618217069'):
    indel_all = pd.DataFrame()
    for sample1 in os.listdir('%s%s/Analysis/INDEL'%(path, s_id)):
        indel_df = pd.read_csv('%s%s/Analysis/INDEL/%s'%(path, s_id, sample1), delimiter='\t')
        indel_df['ID'] = sample1.split('.')[0][6:]
        indel_all = pd.concat([indel_all, indel_df], ignore_index=True, sort=False)
    #indel_all = indel_all.drop(['FORMAT', 'SAMPLE', 'F1R2', 'F2R1', 'FAD', 'SB'], axis=1)
    indel_all.to_html('users_file/%s/Analysis/Results/INDEL.html'%s_id)

    return indel_all

def write_json_cgview(path='users_file/', s_id = 'Test_name_1618217069', baseline_cgv = 1000):
    # For json plots
    
    my_json = json.load(open('static/JS_library/NZ_mito.json'))
    plot_list = []
    track_list = my_json['cgview']['tracks']
    features_list = my_json['cgview']['features']
    path_depths = 'users_file/%s/Analysis/depths/'%s_id
    vcf_list = pd.read_csv('users_file/%s/Analysis/Results/result1.txt'%s_id, delimiter='\t')

    for fastq_file in os.listdir('users_file/%s/fastq/'%(s_id)):
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
    with open('users_file/%s/Analysis/my_cgview.js'%s_id, 'w') as outfile:
        outfile.write(js_str)

    return

def write_json_cgview_indel(path='users_file/', s_id = 'Test_name_1618217069', baseline_cgv = 1000):
    # For json plots
    
    my_json = json.load(open('static/JS_library/NZ_mito.json'))
    plot_list = []
    track_list = my_json['cgview']['tracks']
    features_list = my_json['cgview']['features']
    path_depths = 'users_file/%s/Analysis/depths/'%s_id

    for fastq_file in os.listdir('users_file/%s/fastq/'%(s_id)):
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
        indel_list = read_indel(path='users_file/', s_id = s_id)
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
    with open('users_file/%s/Analysis/my_cgview1.js'%s_id, 'w') as outfile:
        outfile.write(js_str)

    return


from Bio import SeqIO

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

def run_haplogrep3(s_id, vcf_file, tree = 'phylotree-rcrs@17.2', outfile = "haplo.txt"):
    os.system("./haplogrep3/haplogrep3 classify --in %s --tree %s --out %s --extend-report"%(vcf_file, tree, outfile))
    df = pd.read_csv(outfile, delimiter='\t')
    df.to_html('users_file/%s/Analysis/Results/haplogrep.html'%s_id)
    return

