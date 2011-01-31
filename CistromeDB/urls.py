from django.conf.urls.defaults import *

#from django.contrib import admin
from django.contrib.auth.decorators import login_required

from Cistrome.CistromeDB import datacollection
from datacollection.models import Series,Sample

from settings import MEDIA_ROOT

#admin.autodiscover()

urlpatterns = patterns('',
    # Example:
    (r'^accounts/login/$', 'django.contrib.auth.views.login'),
    (r'^accounts/logout/$', 'django.contrib.auth.views.logout'),
    (r'^accounts/profile/$', 'Cistrome.CistromeDB.datacollection.view.userprofile'),
    (r'^acounts/register/$', 'Cistrome.CistromeDB.datacollection.view.register'),
    (r'^CistromeDB/files/(?P<path>.*)$', 'django.views.static.serve', {'document_root': MEDIA_ROOT}),

    (r'^CistromeDB/', include('Cistrome.CistromeDB.datacollection.urls')),

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs' 
    # to INSTALLED_APPS to enable admin documentation:
    #(r'^admin/doc/', include('django.contrib.admindocs.urls')),
    # Uncomment the next line to enable the admin:
    #(r'^$',admin.site.root),
    #(r'^admin/(.*)', admin.site.root),
)
