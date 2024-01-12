from django.shortcuts import render
from django.http import HttpResponse, JsonResponse
from .models import InputForm
from .apps import *
import os
import sys
import mimetypes
from time import time
import hashlib

def index(request):
    os.chdir(os.path.dirname(__file__))
    print(request)
    print(request.FILES)
    if request.method == 'POST':
        form = InputForm(request.POST, request.FILES)
        if form.is_valid():
            if not request.FILES =={}:
                form = form.save(commit=False)
                t1 = int(time())
                t_hash = str(int(time()*100-111223311111)).encode('ASCII')
                m = hashlib.sha256()
                m.update(t_hash)

                session_id = '%s_%s' %(form.name, m.hexdigest())
                session_info = {}
                os.mkdir('users_file/%s' %session_id)
                session_info["organism"] = form.reference_genes
                session_info["analyse_protocol"] = form.analyse_protocol
                session_info["seq_data"]=form.seq_data
                session_info["name"]=form.name
                session_info["email_address"]=form.email_address
                session_info["baseline_cgv"]=form.seq_depth_baseline
                os.mkdir('users_file/%s/fastq'%session_id)
                os.mkdir('users_file/%s/BAMs'%session_id)
                os.mkdir('users_file/%s/Results'%session_id)
                os.mkdir('users_file/%s/Results/QC'%session_id)
                os.mkdir('users_file/%s/Analysis/'%session_id)
                os.mkdir('users_file/%s/Analysis/Results'%session_id)
                sys.stdout = open('users_file/%s/Analysis/Results/log.txt' %session_id, "w")
                handle_uploaded_file(request.FILES['upfile_fastq'], s_id=session_id)
                send_submission(submission_name=form.name, recipient_email=form.email_address)
                print('Upload time %i seconds'%(time()-t1))
                samples = manage_fastq_list(session_id)
                create_yaml(s_id=session_id, samples=samples, info=session_info)
                write_rscript(s_id=session_id)
                return JsonResponse({'progress': 100})
    else:
        form = InputForm()

    return render(request, 'input_mitopore.html', locals())

def download_file():
    print(os.getcwd())
    fl_path = 'mitopore/static/test_result/test.zip'
    filename = 'downloaded_file_name.extension'

    fl = open(fl_path, 'r')
    mime_type, _ = mimetypes.guess_type(fl_path)
    response = HttpResponse(fl, content_type=mime_type)
    
    response['Content-Disposition'] = "attachment; filename=%s" % filename
    return response
