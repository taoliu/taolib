from CistromeDB.datacollection.models import Series, Sample
from django.contrib import admin


class SeriesAdmin(admin.ModelAdmin):
    list_display = ('short_paper','lab','status')
    list_filter = ('status',)
    date_hierarchy = 'pub_date'

admin.site.register(Series,SeriesAdmin)
admin.site.register(Sample)
