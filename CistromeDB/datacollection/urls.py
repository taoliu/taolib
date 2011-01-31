from django.conf.urls.defaults import *
from django.views.generic import list_detail
from django.views.static import serve
from models import Series, Sample
from process import sample_callpeak,sample_ceas
from view import *

series_info = {
    "queryset" : Series.objects.order_by('paper'),
    "paginate_by" : 25,    
    "template_object_name" : "series",
}

samples_info = {
    "queryset" : Sample.objects.all(),
    "paginate_by" : 25,
    "template_object_name" : "samples",
}

urlpatterns = patterns(
    '',
    (r'^$',summary),
    (r'^series/$',list_detail.object_list, series_info),
    (r'^sample/$',list_detail.object_list, samples_info),
    (r'^series/(?P<id>\d+)/sample/$',sample_by_series),
    (r'^series/(?P<series_id>\d+)/$',series_detail),
    (r'^series/(?P<series_id>\d+)/newsample/$',series_add_sample),
    (r'^sample/(?P<sample_id>\d+)/$',sample_detail),
    (r'^sample/species/(?P<species_name>.*?)/$',sample_by_species),
    (r'^sample/cell/(?P<cell_name>.*?)/$',sample_by_cell),
    (r'^sample/exp_type/(?P<exp_type>.*?)/$',sample_by_exptype),
    (r'^sample/platform/(?P<platform>.*?)/$',sample_by_platform),
    (r'^sample/factor/(?P<factor_name>.*?)/$',sample_by_factor),
    (r'^sample/search/$',sample_search),
    (r'^series/search/$',series_search),
    (r'^sample/new/$',sample_new),
    (r'^series/new/$',series_new),
    (r'^sample/(?P<sample_id>\d+)/update/$',sample_update),
    (r'^sample/(?P<sample_id>\d+)/delete/$',sample_delete),
    (r'^sample/(?P<sample_id>\d+)/process/$',sample_process),        
    (r'^series/(?P<series_id>\d+)/update/$',series_update),
    (r'^series/(?P<series_id>\d+)/delete/$',series_delete),
    (r'^comments/', include('django.contrib.comments.urls')),
    (r'^download/(?P<path>.*)$','django.views.static.serve',{'document_root':'/'}), # for development site only!!!
    (r'^sample/(?P<sample_id>\d+)/callpeak/$',sample_callpeak),
    (r'^sample/(?P<sample_id>\d+)/ceas/$',sample_ceas),
    (r'^sample/(?P<sample_id>\d+)/summarycallpeak/$',summary_callpeak),
    (r'^sample/(?P<sample_id>\d+)/summaryceas/$',summary_ceas),
    (r'^sample/(?P<sample_id>\d+)/upload/raw/$',raw_upload),
    (r'^sample/(?P<sample_id>\d+)/download/raw/$',raw_download),    
    (r'^sample/(?P<sample_id>\d+)/upload/peak/$',peak_upload),
    (r'^sample/(?P<sample_id>\d+)/download/peak/$',peak_download),    
    (r'^sample/(?P<sample_id>\d+)/upload/wiggle/$',wiggle_upload),
    (r'^sample/(?P<sample_id>\d+)/download/wiggle/$',wiggle_download),    
)

