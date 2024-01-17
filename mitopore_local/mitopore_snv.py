# from views import *
from apps import manage_fastq_list, create_yaml, run_minimap2, write_json_cgview, write_rscript, run_mutserver, convert_fastq, run_mutect2, run_haplogrep3, write_json_cgview_indel
import os, sys
import time
import shutil
import yaml
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

    #Variance calling using mutserve
    
    # if session_info['analyse_protocol']=='snp_indel':
    #     write_json_cgview_indel(s_id=session_id, baseline_cgv= session_info['baseline_cgv'])
    #     shutil.copyfile('users_file/%s/Analysis/my_cgview1.js'%session_id, 'users_file/%s/JS_library/my_cgview1.js' %session_id)
    #     os.unlink('users_file/%s/Analysis/my_cgview1.js'%session_id)
    #     shutil.copyfile('templates/report.html', 'users_file/%s/report.html'%session_id)
    # else:
    #     shutil.copyfile('templates/report_snv.html', 'users_file/%s/report.html'%session_id)
    # write_json_cgview(s_id=session_id, baseline_cgv= session_info['baseline_cgv'])
    # shutil.copyfile('users_file/%s/Analysis/my_cgview.js'%session_id, 'users_file/%s/JS_library/my_cgview.js' %session_id)
    
    # os.unlink('users_file/%s/Analysis/my_cgview.js'%session_id)
    # os.unlink('users_file/%s/Rplots.pdf'%session_id)

    # if session_info['analyse_protocol']=='snp':
    #     shutil.rmtree('users_file/%s/fastq'%session_id)

    # shutil.make_archive('static/results/%s'%session_id, 'zip', 'users_file/%s/' %session_id)

    # #!TODO: do this in product 

    # link = 'https://www.mitopore.de/static/results/%s.zip'%session_id
    # #send email is not implemented in local mode
    # send_result(str(session_info["name"]), link, session_info["email_address"])
    # try:
    #     shutil.rmtree('users_file/%s'%session_id)
    # except OSError as e:
    #     print ("Error: %s - %s." % (e.filename, e.strerror))
    # return

if __name__ == "__main__":
    main()

# while True:
#     # remove old files
#     directory_path = '/home/ag-rossi/projects/mitopore/mitopore/static/results/'
#     directory_path1 = '/home/ag-rossi/projects/mitopore/mitopore/static/fail/'
#     run_list = os.listdir('users_file')
#     for id1 in run_list:
#         path1= 'users_file/%s/Analysis/submitted.txt'%id1
#         while os.path.isfile(path1):
#             time.sleep(1)
#             #job_execute(session_id=id1)
#             try:
#                 job_execute(session_id=id1)
#             except Exception:
#                 print(Exception)
#                 with open('users_file/%s/config.yaml'%id1) as file:
#                     session_info = yaml.load(file, Loader=yaml.FullLoader)
#                     send_error(submission_name=session_info['name'], recipient_email=session_info['email_address'], err_msg=err_msg)
#                 copy_tree('users_file/%s/'%id1, 'fail/%s'%id1)
#                 shutil.rmtree('users_file/%s'%id1)
#                 pass
#             time.sleep(60)

