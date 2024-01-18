# from views import *
from apps import manage_fastq_list, create_yaml, run_minimap2, write_json_cgview, write_rscript, run_mutserver, convert_fastq, run_mutect2, run_haplogrep3, write_json_cgview_indel
import os, sys
import shutil
from distutils.dir_util import copy_tree

def main():
    path1 = sys.argv[1]
    if path1[-1] == '/':
        path1 = path1[:-1]
    samples = manage_fastq_list(path1)
    create_yaml(path=path1, samples=samples)
    run_minimap2(path=path1)
    run_mutserver(path=path1)
    path2 = os.path.abspath(path1)
    write_rscript(path=path2)
    os.system('R < %s/Analysis/RNA.R --no-save'%path1)
    run_haplogrep3(path=path1)
    write_json_cgview(path=path1)
    copy_tree('static/JS_library/', '%s/JS_library/'%path1)
    shutil.copyfile('%s/Analysis/my_cgview.js'%path1, '%s/JS_library/my_cgview.js' %path1)
    shutil.copyfile('templates/report_snv.html', '%s/report.html'%path1)
    os.unlink('%s/Analysis/my_cgview.js'%path1)
    os.unlink('%s/Rplots.pdf'%path1)

if __name__ == "__main__":
    main()

