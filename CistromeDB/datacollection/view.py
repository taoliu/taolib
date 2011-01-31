from __init__ import *

# modules from python
from shutil import move,rmtree,copyfile
import os
from datetime import datetime
from subprocess import Popen,PIPE

# modules from django
from django.db.models import Q
from django.http import HttpResponse,Http404,HttpResponseRedirect
from django.views.generic import list_detail,create_update
from django.contrib.auth.models import User
from django.shortcuts import get_object_or_404,render_to_response
from django.contrib.auth.forms import UserCreationForm
from django.core.paginator import Paginator, InvalidPage, EmptyPage
from django.core.urlresolvers import reverse
from django.core.files import File
from django.core.files.storage import default_storage

# modules in cistromedb
from models import Sample, Series, Process
import misc


def sample_by_series(request, id):
    # Look up the series id (and raise a 404 if it can't be found).
    series = get_object_or_404(Series,id=id)

    # Use the object_list view for the heavy lifting.
    return list_detail.object_list(
        request,
        queryset = Sample.objects.filter(series=series),
        template_name = "datacollection/sample_list.html",
        template_object_name = "samples",
        paginate_by = 25,
        )

def sample_by_species(request, species_name):
    # Look up the series id (and raise a 404 if it can't be found).

    # Use the object_list view for the heavy lifting.
    return list_detail.object_list(
        request,
        queryset = Sample.objects.filter(species=species_name),
        template_name = "datacollection/sample_species_list.html",
        template_object_name = "samples",
        extra_context={"species":species_name},
        paginate_by = 25,
        )

def sample_by_cell(request, cell_name):
    # Look up the series id (and raise a 404 if it can't be found).

    # Use the object_list view for the heavy lifting.
    return list_detail.object_list(
        request,
        queryset = Sample.objects.filter(cell=cell_name),
        template_name = "datacollection/sample_cell_list.html",
        template_object_name = "samples",
        paginate_by = 25,
        extra_context={"cell":cell_name},
        )

def sample_by_exptype(request, exp_type):
    # Look up the series id (and raise a 404 if it can't be found).

    # Use the object_list view for the heavy lifting.
    return list_detail.object_list(
        request,
        queryset = Sample.objects.filter(exp_type=exp_type),
        template_name = "datacollection/sample_exptype_list.html",
        template_object_name = "samples",
        paginate_by = 25,
        extra_context={"exp_type":exp_type}
        )

def sample_by_platform(request, platform):
    # Look up the series id (and raise a 404 if it can't be found).

    # Use the object_list view for the heavy lifting.
    return list_detail.object_list(
        request,
        queryset = Sample.objects.filter(platform=platform),
        template_name = "datacollection/sample_platform_list.html",
        template_object_name = "samples",
        paginate_by = 25,
        extra_context={"platform":platform}
        )

def sample_by_factor(request, factor_name):
    # Look up the series id (and raise a 404 if it can't be found).

    # Use the object_list view for the heavy lifting.
    return list_detail.object_list(
        request,
        queryset = Sample.objects.filter(factor=factor_name),
        template_name = "datacollection/sample_factor_list.html",
        template_object_name = "samples",
        paginate_by = 25,
        extra_context={"factor":factor_name},
        )

def sample_search(request):
    query = request.GET.get('q','')
    if query:
        qset = (
            Q(cell__icontains=query)|
            Q(factor__icontains=query)|
            Q(species__icontains=query)|
            Q(cell__icontains=query)|
            Q(series__paper__icontains=query)|
            Q(series__lab__icontains=query)
            )
        results = Sample.objects.filter(qset).distinct()
    else:
        results = []

    paginator = Paginator(results, 25)
    try:
        page = int(request.GET.get('page', '1'))
    except ValueError:
        page = 1

    try:
        page_obj = paginator.page(page)
    except (EmptyPage, InvalidPage):
        page_obj = paginator.page(paginator.num_pages)
        
    return list_detail.object_list(
        request,
        queryset = page_obj.object_list,
        template_name = "datacollection/sample_search_list.html",
        template_object_name = "samples",
        extra_context={"query":query,
                       "page_obj":page_obj,
                       "paginator":paginator},
        )

def series_search(request):
    query = request.GET.get('q','')

    if query:
        qset = (
            Q(paper__icontains=query)|
            Q(lab__icontains=query)
            )
        results = Series.objects.filter(qset).distinct()
    else:
        results = []

    paginator = Paginator(results, 25)
    try:
        page = int(request.GET.get('page', '1'))
    except ValueError:
        page = 1

    try:
        page_obj = paginator.page(page)
    except (EmptyPage, InvalidPage):
        page_obj = paginator.page(paginator.num_pages)

        
    return list_detail.object_list(
        request,
        queryset = page_obj.object_list,
        template_name = "datacollection/series_search_list.html",
        template_object_name = "series",
        extra_context={"query":query,
                       "page_obj":page_obj,
                       "paginator":paginator},
        )

def sample_detail(request, sample_id):
    sample = get_object_or_404(Sample,id=sample_id)
    links = __check_available_process (sample)
    
    return list_detail.object_detail(
        request,
        queryset = Sample.objects.all(),
        object_id = sample_id,
        template_name = "datacollection/sample_detail.html",
        template_object_name = "sample",
        extra_context = {"available_link":links}
        )

def __check_available_process (sample):
    # check Peak Calling
    p_callpeak = Process.objects.filter(sample=sample,processtype="P").order_by('-createtime')
    t = []
    if not p_callpeak:
        return [('/CistromeDB/sample/%s/callpeak/' % (str(sample.id)),'Call Peaks'),]
    else:
        p_callpeak = p_callpeak[0]
        p_callpeak.check_status()
        if p_callpeak.status == 'R':
            sample.status='P'
            sample.save()
            return [("",'Call Peaks (running, last check: %s)' % p_callpeak.lastchecktime),]
        elif p_callpeak.status == 'F':
            # move data
            sample.status='C'
            sample.save()
            wd = os.path.join(default_storage.location,'data','raw',str(sample.id))
            peakfile = os.path.join(wd,str(sample.id)+"_peaks.bed")
            wigglefiledir = os.path.join(wd,str(sample.id)+"_MACS_wiggle")            
            if os.path.exists( peakfile ):
                peakdir = os.path.join(default_storage.location,'data','peak',str(sample.id))
                if os.path.exists(peakdir):
                    rmtree(peakdir)
                os.mkdir(peakdir)
                targetfile = os.path.join(peakdir,os.path.basename(peakfile))
                move(peakfile,targetfile )
                sample.peak_file = os.path.join('data','peak',str(sample.id),os.path.basename(peakfile))
                sample.save()
            if os.path.exists( wigglefiledir ):
                wiggledir = os.path.join(default_storage.location,'data','wiggle',str(sample.id))
                if os.path.exists( wiggledir ):
                    rmtree(wiggledir)
                os.mkdir( wiggledir)
                targetfile = os.path.join(wiggledir,os.path.basename(wigglefiledir)+".wig")
                misc.combine(wigglefiledir,targetfile,suffix="gz")
                sample.wiggle_file = os.path.join('data','wiggle',str(sample.id),os.path.basename(wigglefiledir)+".wig")
                sample.save()
                rmtree(wigglefiledir)
            t.append( ("/CistromeDB/sample/%s/summarycallpeak/" % (str(sample.id)),'Call Peaks (finished, %s)' % p_callpeak.finishtime) )
    # check CEAS
    p_ceas = Process.objects.filter(sample=sample,processtype="C").order_by('-createtime')
    if not p_ceas:
        t.append( ('/CistromeDB/sample/%s/ceas/' % (str(sample.id)),'CEAS') )
    else:
        p_ceas = p_ceas[0]
        p_ceas.check_status()
        if p_ceas.status == 'R':
            t.append( ("",'CEAS (running, last check: %s)' % p_ceas.lastchecktime) )
        elif p_ceas.status == 'F':
            t.append( ("/CistromeDB/sample/%s/summaryceas/" % (str(sample.id)),'CEAS (finished, %s)' % p_ceas.finishtime) )
    return t

def summary_ceas (request, sample_id):
    sample = Sample.objects.get(id=sample_id)
    p_ceas = Process.objects.filter(sample=sample,processtype="C").order_by('-createtime')[0]
    if p_ceas.status == 'R':
        return
    else:
        rscript = os.path.join(os.path.dirname(sample.peak_file.path),str(sample_id)+".R")
        os.chdir(os.path.dirname(sample.peak_file.path))
        Popen(["/usr/local/bin/R","--vanilla"],stdin=open(rscript,'r'))
        pdfimg  = os.path.join(os.path.dirname(sample.peak_file.path),str(sample_id)+".pdf")
        pdfdata = open(pdfimg,'r').read()
        response = HttpResponse(pdfdata, mimetype='application/pdf')
        response['Content-Disposition'] = 'attachment; filename=%s' % (str(sample_id)+".pdf")
        return response

def summary_callpeak (request, sample_id):
    sample = Sample.objects.get(id=sample_id)
    p_callpeak = Process.objects.filter(sample=sample,processtype="P").order_by('-createtime')[0]
    if p_callpeak.status == 'R':
        return
    else:
        xls = os.path.join(os.path.dirname(sample.raw_file.path),str(sample_id)+"_peaks.xls")
        xlsdata = open(xls,'r').read()
        response = HttpResponse(xlsdata, mimetype='application/vnd.ms-excel')
        response['Content-Disposition'] = 'attachment; filename=%s' % (str(sample_id)+"_peaks.xls")
        return response


def series_detail(request, series_id):
    series = get_object_or_404(Series,id=series_id)
    samples = Sample.objects.filter(series=series)
    return list_detail.object_detail(
        request,
        queryset = Series.objects.all(),
        object_id = series_id,
        template_name = "datacollection/series_detail.html",
        template_object_name = "series",
        extra_context = {"samples_list":samples,"num_of_samples":len(samples),"num_gt_1":(len(samples)>1)},
        )


def sample_new(request):
    return create_update.create_object(
        request,
        model=Sample,
        login_required=True,
        post_save_redirect="/CistromeDB/sample/%(id)d/upload/raw",
        template_name = "datacollection/new_sample.html",
        extra_context={"user":request.user},
        )

def series_add_sample(request,series_id):
    return create_update.create_object(
        request,
        model=Sample,
        login_required=True,
        post_save_redirect="/CistromeDB/sample/%(id)d/upload/raw",
        template_name = "datacollection/series_new_sample.html",
        extra_context={"user":request.user,"series":Series.objects.get(id=series_id)},
        )

def sample_update(request,sample_id):
    return create_update.update_object(
        request,
        model=Sample,
        login_required=True,
        object_id = sample_id,
        post_save_redirect="/CistromeDB/sample/%(id)d/",
        template_name = "datacollection/update_sample.html",
        extra_context={"user":request.user},
        )

def raw_upload(request,sample_id):
    if request.method == 'POST':
        save_uploaded_file (request.FILES['raw_file'],sample_id,'raw')
        return sample_detail(request,sample_id)
    else:
        return create_update.update_object(
            request,
            model=Sample,
            login_required=True,
            object_id = sample_id,
            post_save_redirect=HttpResponseRedirect(reverse('Cistrome.CistromeDB.datacollection.view.sample_detail', args=(sample_id,))),
            template_name = "datacollection/raw_upload.html",
            extra_context={"user":request.user,},
            )

def save_uploaded_file ( f, sample_id, datatype ):
    sample = Sample.objects.get(id=sample_id)
    if os.path.exists(os.path.join(default_storage.location,'data',datatype,sample_id)):
        rmtree(os.path.join(default_storage.location,'data',datatype,sample_id))
    os.mkdir(os.path.join(default_storage.location,'data',datatype,sample_id))
    destpath = os.path.join(default_storage.location,'data',datatype,sample_id,f.name)
    destname = os.path.join('data',datatype,sample_id,f.name)
    
    if f.multiple_chunks():
        # big file
        move(f.temporary_file_path(),destpath)
    else:
        # small file
        ff = open(destpath,'wb')
        ff.write(f.read())
        ff.close
    sample.raw_file = destname
    sample.save()
    return
    

def peak_upload(request,sample_id):
    return create_update.update_object(
        request,
        model=Sample,
        login_required=True,
        object_id = sample_id,
        post_save_redirect="/CistromeDB/sample/%(id)d/",
        template_name = "datacollection/peak_upload.html",
        extra_context={"user":request.user},
        )

def wiggle_upload(request,sample_id):
    return create_update.update_object(
        request,
        model=Sample,
        login_required=True,
        object_id = sample_id,
        post_save_redirect="/CistromeDB/sample/%(id)d/",
        template_name = "datacollection/wiggle_upload.html",
        extra_context={"user":request.user},
        )

def raw_download(request,sample_id):
    return

def peak_download(request,sample_id):
    return

def wiggle_download(request,sample_id):
    return
    
def sample_delete(request,sample_id):
    return create_update.delete_object(
        request,
        model=Sample,
        login_required=True,
        object_id=sample_id,
        post_delete_redirect="/CistromeDB/sample/",
        template_name = "datacollection/delete_sample_confirmation.html",
        extra_context={"user":request.user},
        )

def series_new(request):
    return create_update.create_object(
        request,
        model=Series,
        login_required=True,        
        template_name = "datacollection/new_series.html",
        post_save_redirect="/CistromeDB/series/%(id)d/",
        extra_context={"user":request.user},        
        )

def series_update(request,series_id):
    return create_update.update_object(
        request,
        model=Series,
        login_required=True,
        object_id = series_id,
        post_save_redirect="/CistromeDB/series/%(id)d/",
        template_name = "datacollection/update_series.html",
        extra_context={"user":request.user},        
        )

def series_delete(request,series_id):
    return create_update.delete_object(
        request,
        model=Series,
        login_required=True,
        object_id=series_id,
        post_delete_redirect="/CistromeDB/series/",
        template_name = "datacollection/delete_series_confirmation.html",
        extra_context={"user":request.user},        
        )


def userprofile(request):
    return render_to_response("datacollection/userprofile.html",{
        "user"   : request.user,
        })

def register(request):
    return render_to_response("datacollection/userprofile.html",{
        "user"   : request.user,
        })

def sample_process (request,sample_id):
    sample = get_object_or_404(Sample,id=sample_id)
    (s_raw_file,s_peak_file,s_wiggle_file) = (sample.raw_file,sample.peak_file,sample.wiggle_file)
    dir1 = os.path.join(REPO_DIR,sample_id)
    if os.path.exists( dir1 ):
        rmtree( dir1 )
    
    os.mkdir( dir1 )
    os.mkdir( os.path.join(dir1,'raw') )
    os.mkdir( os.path.join(dir1,'peak') )
    os.mkdir( os.path.join(dir1,'wiggle') )        

    if os.path.isfile(sample.raw_file) and os.path.isfile(sample.peak_file) and os.path.isfile(sample.wiggle_file):
        incomplete = False
    else:
        incomplete = True

    if os.path.isfile(sample.raw_file):
        t_raw_file = os.path.join(dir1,'raw',os.path.basename(s_raw_file))
        if os.path.basename(s_raw_file) == 'NA':
            copyfile (s_raw_file,t_raw_file)
        else:
            move (s_raw_file,t_raw_file)
        sample.raw_file = t_raw_file
        
    if os.path.isfile(sample.peak_file):
        t_peak_file = os.path.join(dir1,'peak',os.path.basename(s_peak_file))        
        if os.path.basename(s_peak_file) == 'NA':
            copyfile (s_peak_file,t_peak_file)
        else:
            move (s_peak_file,t_peak_file)
        sample.peak_file = t_peak_file
        
    if os.path.isfile(sample.wiggle_file):
        t_wiggle_file = os.path.join(dir1,'wiggle',os.path.basename(s_wiggle_file))
        if os.path.basename(s_wiggle_file) == 'NA':
            copyfile (s_wiggle_file,t_wiggle_file)
        else:
            move (s_wiggle_file,t_wiggle_file)
        sample.wiggle_file = t_wiggle_file
    
    if incomplete:
        sample.status = 'R'
        sample.save()
        return render_to_response("datacollection/sample_detail_process.html",{
            "sample" : sample,
            })
    else:
        sample.status = 'P'
        sample.save()
        return render_to_response("datacollection/sample_detail_process.html",{
            "sample" : sample,
            })

def summary(request):
    series = Series.objects.all()
    number_of_series = len(series)
    samples = Sample.objects.all()
    number_of_sample = len(samples)    
    return render_to_response("datacollection/summary.html",{
        "nseries"   : number_of_series,
        "series"    : series[:5],
        "nsample"   : number_of_sample,
        "samples"    : samples[:5],       
        "now"       : datetime.now(),
        })
