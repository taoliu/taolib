from __init__ import *
from django.db import models
from django.utils.encoding import smart_unicode
from django.contrib.auth.models import User
from datetime import datetime
from subprocess import Popen,PIPE

# Create your models here.
SPECIES_CHOICES = (
    ('human', 'Homo sapiens'),
    ('fruit fly'  , 'Drosophila melanogaster'),
    ('nematode' , 'Caenorhabditis elegans'),
    ('house mouse', 'Mus musculus'),
    ('black rat'  , 'Rattus rattus')
    )

EXP_TYPE_CHOICES = (
    ('array','ChIP-chip'),
    ('seq','ChIP-seq'),
    )

STATUS_CHOICES = (
    ('C','Completed'),
    ('P','In Processing'),
    ('R','Request data'),
    ('O','Obsolete'),
    )

CHIP_TYPE_CHOICES = (
    ('promoter','Promoter Array'),
    ('wg','Whole Genome Array'),
    )

CHIP_MANUFACTURE = (
    ('affy','Affymetrix'),
    ('nimble','NimbleGen'),
    ('agilent','Agilent'),
    )

SEQUENCING_TYPE_CHOICES = (
    ('solexa',''),
    ('solid',''),
    ('helicos'),
    ('mpss',''),
    ('454',''),
    ('sanger',''),
    ('other',''),
    )

PROCESS_STATUS_CHOICES = (
    ('F','Finished'),
    ('R','Running'),
    ('E','Error'),
    )

PROCESS_TYPE_CHOICES = (
    ('P','Call Peaks'),
    ('C','CEAS'),
    ('S','SEQPOS'),
    )

class Series(models.Model):
    """For serie of experiments linked to a single publication.
    
    """
    paper = models.CharField(max_length=300,
                             help_text="Title of paper or simple description")
    lab   = models.CharField(max_length=300,default="NA",
                             help_text="Name of the lab/corresponding author")
    pubmed = models.IntegerField(default=0,
                                 help_text="Pubmed id, default: 0")
    geo    = models.IntegerField(default=0,
                                 help_text="GEO series id, for GSE10023, input 10023, default: 0")
    pub_date = models.DateField('Publication Date',
                                help_text="date of publication. Format:YYYY-MM-DD")
    status = models.CharField(max_length=50,
                              choices=STATUS_CHOICES)
    def __unicode__ (self):
        return self.short_paper()

    def short_paper (self):
        v = smart_unicode(self.paper)
        if len(v) > 50:
            s = v[:47]+u'...'
        else:
            s = v
        return s
    
    short_paper.short_description =  "Paper"

class Sample(models.Model):
    """For a sample in a series.
    
    """
    series = models.ForeignKey(Series)
    factor = models.CharField(max_length=200,help_text="name of transcription factor/histone modification")
    species = models.CharField(max_length=200,choices=SPECIES_CHOICES)
    cell  = models.CharField(max_length=200,default="NA")
    #cell = models.ForeignKey(Cell)
    exp_type = models.CharField(max_length=200,choices=EXP_TYPE_CHOICES,
                                help_text="type of experiment")
    platform = models.CharField(max_length=300,default="NA")
    status = models.CharField(max_length=50,
                              choices=STATUS_CHOICES)
    raw_file = models.FileField(upload_to="data/raw",help_text="Select the raw file you want to upload",blank=True)
    peak_file = models.FileField(upload_to="data/peak",help_text="Select the peak file you want to upload",blank=True)
    wiggle_file = models.FileField(upload_to="data/wiggle",help_text="Select the wiggle file you want to upload",blank=True)
    
    def __unicode__ (self):
        return self.description()

    def description (self):
        f = smart_unicode(self.factor)
        s = smart_unicode(self.species)
        c = smart_unicode(self.cell)
        p = smart_unicode(self.platform)
        return f+" IN "+s+" "+c+" CELL" + " BY " + p

class Cell(models.Model):
    """Cell line information.
    
    """
    name = models.CharField(max_length=100)
    species = models.CharField(max_length=200,choices=SPECIES_CHOICES)
    description = models.TextField()

    def __unicode__ (self):
        return self.name()

    def short_description (self):
        v = smart_unicode(self.description)
        if len(v) > 50:
            s = v[:47]+u'...'
        else:
            s = v
        return s
    
class DBUserProfile(models.Model):
    user = models.ForeignKey(User, unique=True)
    
class TFSummary(models.Model):
    sample = models.ForeignKey(Sample)
    ceas_file = models.CharField(max_length=300)
    

class Process(models.Model):
    """Record process information.
    
    """
    user = models.CharField(max_length=255)
    sample = models.ForeignKey(Sample)
    qsubpid = models.CharField(max_length=255)
    createtime = models.DateTimeField(auto_now_add=True)
    lastchecktime = models.DateTimeField(auto_now=True)
    finishtime = models.DateTimeField(null=True)
    status = models.CharField(max_length=50,
                              choices=PROCESS_STATUS_CHOICES)
    processtype = models.CharField(max_length=50,
                                   choices=PROCESS_TYPE_CHOICES)

    def check_status (self):
        statusinfo = Popen([QSTAT, "-r", str(self.qsubpid)], stdout=PIPE).communicate()[0]
        if not statusinfo:
            self.status = 'F'
            if not self.finishtime:
                self.finishtime = datetime.now()
            self.save()
        else:
            self.status = 'R'
            self.save()

