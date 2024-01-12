# from views import *
from apps import run_minimap2, send_result, send_error, write_json_cgview, write_rscript, run_mutserver, convert_fastq, run_mutect2, run_haplogrep3, write_json_cgview_indel
import os
import time
import shutil
import yaml
from distutils.dir_util import copy_tree
# from django.core import mail
# from django.core.mail.backends.smtp import EmailBackend
import datetime
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'settings')
err_msg = ''
def remove_old_files(directory):
    now = datetime.datetime.now()
    threshold = now - datetime.timedelta(days=7)

    for root, dirs, files in os.walk(directory, topdown=False):
        for file_name in files:
            file_path = os.path.join(root, file_name)
            modification_time = datetime.datetime.fromtimestamp(os.path.getmtime(file_path))
            if modification_time < threshold:
                os.remove(file_path)
                print(f"Removed file: {file_path}")

        for dir_name in dirs:
            dir_path = os.path.join(root, dir_name)
            modification_time = datetime.datetime.fromtimestamp(os.path.getmtime(dir_path))
            if modification_time < threshold:
                os.rmdir(dir_path)
                print(f"Removed directory: {dir_path}")

def job_execute(session_id):
    print('Start running job for task %s'%session_id)
    with open('users_file/%s/config.yaml'%session_id) as file:
        session_info = yaml.load(file, Loader=yaml.FullLoader)

    if session_info['analyse_protocol']=='snp_indel':
        for files in os.listdir('users_file/%s/fastq/'%session_id):
            if not('bam' in files):
                convert_fastq(fastqfile='users_file/%s/fastq/%s'%(session_id, files), outfile='users_file/%s/fastq/%s'%(session_id, files))
                print('convert1')
            else:
                err_msg = 'Indel pipeline use fastq files only'
        run_mutect2(s_id = session_id, organism=session_info["organism"])

    elif not('bam' in os.listdir('users_file/%s/fastq/'%session_id)[0]):   
            t3 = time.time() 
            run_minimap2(s_id=session_id, seq_data=session_info["seq_data"], organism=session_info["organism"])
            t4 = time.time()
            print('Run Minimap time %i seconds' %(t4-t3))
    else:
        if not os.path.exists('users_file/%s/Analysis'%session_id):
            os.mkdir('users_file/%s/Analysis'%session_id)
        path_minimap = 'users_file/%s/Analysis/Minimap'%session_id
        if not os.path.exists(path_minimap):
            os.mkdir(path_minimap)
        # Run R for quality control
        for files in os.listdir('users_file/%s/fastq/'%session_id):
            if 'bam' in files:
                shutil.move('users_file/%s/fastq/%s'%(session_id, files), 'users_file/%s/Analysis/Minimap/'%session_id)   
    run_mutserver(s_id=session_id)

    os.system('R < users_file/%s/Analysis/RNA.R --no-save'%session_id)
    run_haplogrep3(s_id=session_id, vcf_file= 'users_file/%s/Analysis/Results/result1.vcf'%session_id, outfile='users_file/%s/Analysis/Results/haplogrep3.txt'%session_id)
    copy_tree('static/JS_library/', 'users_file/%s/JS_library/'%session_id)

    #Variance calling using mutserve
    
    if session_info['analyse_protocol']=='snp_indel':
        write_json_cgview_indel(s_id=session_id, baseline_cgv= session_info['baseline_cgv'])
        shutil.copyfile('users_file/%s/Analysis/my_cgview1.js'%session_id, 'users_file/%s/JS_library/my_cgview1.js' %session_id)
        os.unlink('users_file/%s/Analysis/my_cgview1.js'%session_id)
        shutil.copyfile('templates/report.html', 'users_file/%s/report.html'%session_id)
    else:
        shutil.copyfile('templates/report_snv.html', 'users_file/%s/report.html'%session_id)
    write_json_cgview(s_id=session_id, baseline_cgv= session_info['baseline_cgv'])
    shutil.copyfile('users_file/%s/Analysis/my_cgview.js'%session_id, 'users_file/%s/JS_library/my_cgview.js' %session_id)
    
    os.unlink('users_file/%s/Analysis/my_cgview.js'%session_id)
    os.unlink('users_file/%s/Rplots.pdf'%session_id)

    if session_info['analyse_protocol']=='snp':
        shutil.rmtree('users_file/%s/fastq'%session_id)

    shutil.make_archive('static/results/%s'%session_id, 'zip', 'users_file/%s/' %session_id)

    #!TODO: do this in product 

    link = 'https://www.mitopore.de/static/results/%s.zip'%session_id
    #send email is not implemented in local mode
    send_result(str(session_info["name"]), link, session_info["email_address"])
    try:
        shutil.rmtree('users_file/%s'%session_id)
    except OSError as e:
        print ("Error: %s - %s." % (e.filename, e.strerror))
    return


while True:
    # remove old files
    directory_path = '/home/ag-rossi/projects/mitopore/mitopore/static/results/'
    directory_path1 = '/home/ag-rossi/projects/mitopore/mitopore/static/fail/'
    remove_old_files(directory_path)
    remove_old_files(directory_path1)
    run_list = os.listdir('users_file')
    for id1 in run_list:
        path1= 'users_file/%s/Analysis/submitted.txt'%id1
        while os.path.isfile(path1):
            time.sleep(1)
            #job_execute(session_id=id1)
            try:
                job_execute(session_id=id1)
            except Exception:
                print(Exception)
                with open('users_file/%s/config.yaml'%id1) as file:
                    session_info = yaml.load(file, Loader=yaml.FullLoader)
                    send_error(submission_name=session_info['name'], recipient_email=session_info['email_address'], err_msg=err_msg)
                copy_tree('users_file/%s/'%id1, 'fail/%s'%id1)
                shutil.rmtree('users_file/%s'%id1)
                pass
            time.sleep(60)

