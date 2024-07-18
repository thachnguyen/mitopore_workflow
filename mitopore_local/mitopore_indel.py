# from views import *
from apps import manage_fastq_list, create_yaml, run_minimap2, write_json_cgview, write_rscript, run_mutserver, convert_fastq, run_mutect2, run_haplogrep3, write_json_cgview_indel
import os, sys
import shutil
from distutils.dir_util import copy_tree

def main():
    os.chdir('/home/ag-rossi/projects/mitopore_workflow/mitopore_local/')
    if len(sys.argv) > 1:
        path1 = sys.argv[1]
        if path1[-1] == '/':
            path1 = path1[:-1]
    else: 
        path1 = '/mitopore_data'
    samples = manage_fastq_list(path1)
    create_yaml(path=path1, samples=samples)

    for files in os.listdir('%s/fastq/'%path1):
        if 'fastq' == files[-5:]:
            if len(sys.argv) > 2:
                if sys.argv[2]!='illumina':
                    if not os.path.exists('%s/fastq_fil/'%path1):
                        os.mkdir('%s/fastq_fil/'%path1)
                    convert_fastq(fastqfile='%s/fastq/%s'%(path1, files), outfile='%s/fastq_fil/%s'%(path1, files))
                    print('converting fastq using Elibq')
                else:
                    print('Ilumina')
                    if os.path.exists('%s/fastq_fil/'%path1):
                        # Remove the destination directory
                        shutil.rmtree('%s/fastq_fil/'%path1)
                    shutil.copytree('%s/fastq/'%path1, '%s/fastq_fil/'%path1)
                    
            else:
                    if not os.path.exists('%s/fastq_fil/'%path1):
                        os.mkdir('%s/fastq_fil/'%path1)
                    convert_fastq(fastqfile='%s/fastq/%s'%(path1, files), outfile='%s/fastq_fil/%s'%(path1, files))
                    print('converting fastq using Elibq')
        else:
            print('Indel pipeline use fastq files only')
    run_mutect2(path=path1)


    run_minimap2(path=path1)
    run_mutserver(path=path1)
    path2 = os.path.abspath(path1)
    write_rscript(path=path2)
    os.system('R < %s/Analysis/RNA.R --no-save'%path1)
    run_haplogrep3(path=path1)
    write_json_cgview(path=path1)
    write_json_cgview_indel(path=path1)

    copy_tree('static/JS_library/', '%s/JS_library/'%path1)
    shutil.copyfile('%s/Analysis/my_cgview.js'%path1, '%s/JS_library/my_cgview.js' %path1)
    shutil.copyfile('%s/Analysis/my_cgview1.js'%path1, '%s/JS_library/my_cgview1.js' %path1)
    shutil.copyfile('templates/report.html', '%s/report.html'%path1)
    os.unlink('%s/Analysis/my_cgview.js'%path1)
    os.unlink('%s/Analysis/my_cgview1.js'%path1)
    os.unlink('%s/Rplots.pdf'%path1)

if __name__ == "__main__":
    main()