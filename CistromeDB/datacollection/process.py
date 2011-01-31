from django.http import Http404,HttpResponseRedirect,HttpResponse
from django.core.files.storage import default_storage
from models import Sample, Series, Process

from __init__ import *

from shutil import move,rmtree,copyfile
import os
from datetime import datetime

from django.contrib.auth.models import User
from django.shortcuts import get_object_or_404,render_to_response

from django import forms

from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse

import tempfile
from subprocess import Popen,PIPE

DATAPATH = default_storage.location

EFFECTIVE={
    'human':2700000000,
    'fruit fly':100000000,
    'nematode':100000000,
    'house mouse':2000000000,
    'black rat':2000000000,
    }

GENOMEASSEMBLE={
    'human':'hg18',
    'house mouse':'mm8',
    'fruit fly':'dm3',
    'nematode':'ce4',
    }

MACSFORMATCHOICES=(
    ("BED","UCSC BED format"),
    ("ELAND","ELAND result format"),
    ("ELANDMULTI","ELAND multiple alignment format"),
)

class MACSForm (forms.Form):
    format= forms.ChoiceField(choices=MACSFORMATCHOICES,help_text='Format of tag file, "ELAND" or "BED" or "ELANDMULTI" or "ELANDMULTIPET".')
    mfold = forms.IntegerField(initial=32,help_text='Select the regions with MFOLD high-confidence enrichment ratio against background to build model.')
    tsize = forms.IntegerField(initial=25,help_text='Tag size.')
    gsize = forms.IntegerField(required=False,help_text='Effective genome size.')
    bw    = forms.IntegerField(initial=300,help_text='Band width. This value is used while building the shifting model.')
    pvalue= forms.FloatField(initial=1e-5,help_text='Pvalue cutoff for peak detection.')
    nolmbd= forms.BooleanField(required=False,initial=False,help_text='If True, MACS will use fixed background lambda as local lambda for every peak region. Normally, MACS calculates a dynamic local lambda to reflect the local bias due to potential chromatin structure.')

class MA2CForm (forms.Form):
    pvalue= forms.FloatField()

@login_required
def sample_callpeak (request,sample_id):
    """Call peaks from Raw data for a sample.
    
    """
    sample = Sample.objects.get(id=sample_id)
    if sample.exp_type == 'seq':
        if request.method == 'POST':
            form = MACSForm(request.POST)
            if form.is_valid():
                return __run_macs( request,sample, form)
        else:
            form = MACSForm()
        return render_to_response('datacollection/runmacs.html',{
                'form':form,
                })
    elif sample.exp_type == 'chip':
        if request.method == 'POST':
            form = MA2CForm(request.POST)
            if form.is_valid():
                return __run_ma2c( request,sample, form)
        else:
            form = MA2CForm()
        return reder_to_response('datacollection/runma2c.html',{
                'form':form,
                })

def __run_macs ( request,sample, form):
    script_file_name = __create_macs_script(sample,form)
    return __submit_qsub (request,script_file_name,sample,"P")

def __create_macs_script ( sample,form ):
    dirpath = os.path.dirname(sample.raw_file.path)
    scriptname = os.path.basename(tempfile.NamedTemporaryFile('w').name)
    scriptpath = os.path.join(dirpath,'callpeak_'+scriptname+'.sh')
    f = open(scriptpath,'w')
    os.chmod(scriptpath,0666)
    rawfile = os.path.basename(sample.raw_file.path)
    eff = EFFECTIVE[sample.species]
    f.write("""#PBS -S /bin/bash
#PBS -N ccp%s
#PBS -l walltime=1000:00:00,nodes=1:ppn=1
#PBS -k oe
#PBS -o %s
#PBS -e %s
""" % (f.name,dirpath,dirpath))
    f.write("cd %s\n" % dirpath)
    f.write("%s -t %s --name %s --format %s --gsize %d --wig --tsize %d --bw %d --mfold %d --pvalue %e " %
            (MACS,rawfile,str(sample.id),form.cleaned_data['format'],eff,
             form.cleaned_data['tsize'],form.cleaned_data['bw'],
             form.cleaned_data['mfold'],form.cleaned_data['pvalue'],
             ))
    if form.cleaned_data['nolmbd']:
        f.write('--nolambda\n')
    else:
        f.write('\n')
    f.close()
    return f.name

def __create_callpeak_script_ma2c (sample):
    pass

@login_required
def sample_ceas ( request, sample_id ):
    sample = Sample.objects.get(id=sample_id)
    script_file_name = __create_ceas_script(sample)
    return __submit_qsub (script_file_name,sample_id,request,"C")

class CEASForm(forms.Form):
    chip_res   = forms.IntegerField(initial=600,min_value=600,help_text="ChIP annotation resolution. Minimum Value is 600bp.")
    promoter   = forms.IntegerField(initial=3000,max_value=10000,help_text="Promoter size for annotation. Maximum Value is 10000bp.")
    bipromoter = forms.IntegerField(initial=5000,max_value=20000,help_text="Bidirectional-promoter size for annotation. Maximum Value is 20000bp.")
    downstream = forms.IntegerField(initial=3000,max_value=10000,help_text="Downstream size for annotation. Maximum Value is 10000bp.")
    pf_res     = forms.IntegerField(initial=50,min_value=10,help_text="Wig profiling resolution. Minimum Value is 10bp.")
    rel_dist   = forms.IntegerField(initial=3000,help_text="Relative distance to TSS/TTS in wig profiling.")
    metagene   = forms.IntegerField(initial=3000,min_value=10,help_text="Wig profiling resolution. Minimum Value is 10bp.") 

def __create_ceas_script(sample):
    dirpath = os.path.dirname(sample.peak_file.path)
    scriptname = os.path.basename(tempfile.NamedTemporaryFile('w').name)
    scriptpath = os.path.join(dirpath,'ceas_'+scriptname+'.sh')
    f = open(scriptpath,'w')
    os.chmod(scriptpath,0666)

    gtname = os.path.join(GENOMEASSEMBLE[sample.species])
    # ceas.py  -g ce4 -w diPETnomodel.wig -b diPETmodel_peaks.bed --name diPETnomodel --rel-dist=1000 --bg
    f.write("""#PBS -S /bin/bash
#PBS -N cceas%s
#PBS -l walltime=1000:00:00,nodes=1:ppn=1
#PBS -k oe
#PBS -o %s
#PBS -e %s
""" % (f.name,dirpath,dirpath))
    f.write("cd %s\n" % dirpath)
    f.write("%s -g %s --bg -b %s -w %s --name %s --rel-dist=3000\n" %
            (CEAS,os.path.join(GENETABLEDIR,gtname),sample.peak_file.path,
             sample.wiggle_file.path,str(sample.id)))
    f.close()
    return f.name
    
def __submit_qsub( request,script, sample,processtype ):
    pid = Popen ([QSUB,script],stdout=PIPE).communicate()[0].strip()
    # 5490.freesia.dfci.harvard.edu
    process = Process.objects.create(user=request.user.username,sample=sample,qsubpid=pid,status="R",processtype=processtype)
    return HttpResponseRedirect(reverse('Cistrome.CistromeDB.datacollection.view.sample_detail', args=(sample.id,)))
