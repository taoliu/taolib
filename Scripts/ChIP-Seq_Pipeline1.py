#!/usr/bin/env python
# Time-stamp: <2011-03-08 16:33:11 Tao Liu>

"""Script Description: A demo ChIP-seq pipeline script. From reads
mapping to motif analysis. It will do 4 validity checks before running
the pipeline.

Copyright (c) 2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import re
from optparse import OptionParser
import ConfigParser
import string
import logging
from subprocess import call as subpcall
from os.path import join as pjoin

# ------------------------------------
# constants
# ------------------------------------

logfhd = open("log","w")

logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )

error   = logging.critical		# function alias
warn    = logging.warning

def info(a):
    logging.info(a)
    logfhd.write(a+"\n")
    logfhd.flush()

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

def read_config(configFile):
    """Read configuration file and parse it into a dictionary.

    In the dictionary, the key is the section name plus option name like: data.data.treatment_seq_file_path.
    """
    configs = {}
    
    configParser = ConfigParser.ConfigParser()

    if len(configParser.read(configFile)) == 0:
        raise IOError("%s not found!" % configFile)
    
    for sec in configParser.sections():
        secName = string.lower(sec)
        for opt in configParser.options(sec):
            optName = string.lower(opt)
            configs[secName + "." + optName] = string.strip(configParser.get(sec, opt))
    
    return configs

def check_conf_validity(config_dict):
    """A configuration validity test. Return false if any required option is not valid.

    """
    if not config_dict["sample.sample_id"].isdigit():
        error("Config Validity Check: sample id -- %s should be specified and be an integer!" % config_dict["sample.sample_id"])
        return False
    if not config_dict["sample.assembly_name"]:
        error("Config Validity Check: assembly name -- %s should be specified!")
        return False
    else:
        an = config_dict["sample.assembly_name"]
        if not an.startswith("hg") and not an.startswith("mm"):
            error("Config Validity Check: assembly name -- %s should be human or mouse!" % an)
            return False
        elif not an[2:].isdigit():
            error("Config Validity Check: assembly name -- %s should be UCSC dbkey!" % an)
            return False
        elif an.startswith("hg"):
            config_dict["sample.species"] = "hs"
            info("Species is 'human'")            
        elif an.startswith("mm"):
            config_dict["sample.species"] = "mm"
            info("Species is 'mouse'")
            
    rdt = config_dict["data.raw_data_type"]
    if rdt != "seq" and rdt != "alignment": # and rdt != "peakcalls":
        error("Config Validity Check: raw data type -- %s should be either seq, alignment!" % rdt)
        return False

    rdf = config_dict["data.raw_data_format"]
    # check if the raw data format matches raw data type
    if rdt == "seq":
        if ( rdf != "fastq" and rdf != "fasta" and rdf != "sequences"):
            error("Config Validity Check: raw data is seq, then data format -- %s should be either fastq, fasta or sequences!" % rdf)
            return False            
        elif ( config_dict["data.treatment_seq_file_path"] == "" ): #or config_dict["data.control_seq_file_path"] == "" ):
            error("Config Validity Check: raw data is seq however sequencing raw files are missing!" )
            return False
    elif rdt == "alignment":
        if ( rdf != "SAM" and rdf != "BAM" and rdf != "BED"):
            error("Config Validity Check: raw data is alignment, then data format -- %s should be either SAM, BAM or BED!" % rdf)
            return False
        elif ( config_dict["data.treatment_ali_file_path"] == "" ): #or config_dict["data.control_ali_file_path"] == "" ):
            error("Config Validity Check: raw data is alignment however alignment raw files are missing!" )
            return False

    # other required options:
    if not config_dict["bowtie.bowtie_main"]:
        error("Config Validity Check: bowtie main is missing!")
        return False
    if not config_dict["bowtie.bowtie_genome_index_path"]:
        error("Config Validity Check: bowtie genome index is missing!")
        return False
    if not config_dict["samtools.samtools_main"]:
        error("Config Validity Check: samtools main is missing!")
        return False
    if not config_dict["samtools.samtools_chrom_len_path"]:
        error("Config Validity Check: samtools chromosome length file is missing!")
        return False    
    if not config_dict["macs.macs_main"]:
        error("Config Validity Check: macs main is missing!")
        return False
    if not config_dict["ceas.ceas_main"]:
        error("Config Validity Check: ceas main is missing!")
        return False
    if not config_dict["ceas.ceas_genetable_path"]:
        error("Config Validity Check: ceas genetable is missing!")
        return False
    if not config_dict["venn.venn_diagram_main"]:
        error("Config Validity Check: Venn diagram main is missing!")
        return False
    if not config_dict["venn.dhs_bed_path"]:
        error("Config Validity Check: DHS bed file is missing!")
        return False
    if not config_dict["correlation.wig_correlation_main"]:
        error("Config Validity Check: wig correlation main is missing!")
        return False
    avail = ["mean","median","sample"]
    wcm = config_dict["correlation.wig_correlation_method"]
    if avail.count(wcm) == 0:
        error("Config Validity Check: Correlation method '%s' is invalid! It should be either mean, median, or sample!" % wcm)
    if not config_dict["conservation.conserv_plot_main"]:
        error("Config Validity Check: conservation plot main is missing!")
        return False
    if not config_dict["conservation.conserv_plot_phast_path"]:
        error("Config Validity Check: PhastCons data file folder is missing!")
        return False
    if not config_dict["seqpos.seqpos_main"]:
        error("Config Validity Check: seqpos main is missing!")
        return False
    #if not config_dict["seqpos.seqpos_motif_db_path"]:
    #    error("Config Validity Check: seqpos motif db path is missing!")
    #    return False
    sms = config_dict["seqpos.seqpos_motif_db_selection"].split(",")
    avail = ["transfac.xml","pbm.xml", "jaspar.xml", "hpdi.xml", "y1h.xml"]
    for sm in sms:
        if avail.count(sm) == 0:
            error("Config Validity Check: seqpos motif db %s is invalid! It should be either transfac.xml, pbm.xml, jaspar.xml, hpdi.xml, or y1h.xml!" % sm)
            return False
    return True

def check_file_validity(config_dict):
    """A data file validity test. Return false if any file can't be found.

    """
    if config_dict["data.control_seq_file_path"]:
        config_dict["data.has_control"] = True
        info("Data has control.")
    else:
        config_dict["data.has_control"] = False
    if config_dict["data.raw_data_type"] == 'seq':
        # check the raw data file option
        tfiles = config_dict["data.treatment_seq_file_path"].split(",")
        for tfile in tfiles:
            if not os.path.isfile(tfile):
                error("File Validity Check: %s can't be found!" % tfile)
                return False
        if config_dict["data.has_control"]:
            cfiles = config_dict["data.control_seq_file_path"].split(",")
            for cfile in cfiles:
                if not os.path.isfile(cfile):
                    error("File Validity Check: %s can't be found!" % cfile)
                    return False
    elif config_dict["data.raw_data_type"] == 'alignment':
        # check the alignment data file option
        tfiles = config_dict["data.treatment_ali_file_path"].split(",")
        for tfile in tfiles:
            if not os.path.isfile(tfile):
                error("File Validity Check: %s can't be found!" % tfile)
                return False
        if config_dict["data.has_control"]:            
            cfiles = config_dict["data.control_ali_file_path"].split(",")
            for cfile in cfiles:
                if not os.path.isfile(cfile):
                    error("File Validity Check: %s can't be found!" % cfile)
                    return False

    return True

def check_cmd_validity (config_dict):
    """Check if command line can be found.
    
    """
    if not os.path.isfile(config_dict["bowtie.bowtie_main"]):
       error("CMD Validity Check: bowtie main can't be found!")
       return False
    if not os.path.isfile(config_dict["samtools.samtools_main"]):
       error("CMD Validity Check: samtools main can't be found!")
       return False
    if not os.path.isfile(config_dict["macs.macs_main"]):
        error("CMD Validity Check: macs main can't be found!")
        return False
    if not os.path.isfile(config_dict["ceas.ceas_main"]):
        error("CMD Validity Check: ceas main can't be found!")
        return False
    if not os.path.isfile(config_dict["venn.venn_diagram_main"]):
       error("CMD Validity Check: Venn diagram main can't be found!")
       return False
    if not os.path.isfile(config_dict["correlation.wig_correlation_main"]):
        error("CMD Validity Check: wig correlation main can't be found!")
        return False
    if not os.path.isfile(config_dict["conservation.conserv_plot_main"]):
        error("CMD Validity Check: conservation plot main can't be found!")
        return False
    if not os.path.isfile(config_dict["seqpos.seqpos_main"]):
       error("CMD Validity Check: seqpos main can't be found!")
       return False
    return True

def check_lib_validity (config_dict):
    """Check if the library files for command line can be found.
    
    """
    if not os.path.isfile(config_dict["bowtie.bowtie_genome_index_path"]+".1.ebwt"):
        error("Library Files Validity Check: bowtie genome index can't be found!")
        return False
    if not os.path.isfile(config_dict["samtools.samtools_chrom_len_path"]):
        error("Library Files Validity Check: samtools chromosome length file can't be found!")
        return False    
    if not os.path.isfile(config_dict["ceas.ceas_genetable_path"]):
        error("Library Files Validity Check: ceas genetable can't be found!")
        return False
    if not os.path.isfile(config_dict["venn.dhs_bed_path"]):
        error("Library Files Validity Check: DHS bed file can't be found!")
        return False
    if not os.path.isdir(config_dict["conservation.conserv_plot_phast_path"]):
        error("Library Files Validity Check: PhastCons data file folder can't be found!")
        return False
    #sms = config_dict["seqpos.seqpos_motif_db_selection"].split(",")
    #for sm in sms:
    #    if not os.path.isfile(pjoin(config_dict["seqpos.seqpos_motif_db_path"],sm)):
    #        error("Library Files Validity Check: seqpos motif db %s can't be found!" % sm)
    #        return False
    return True

def prepare_output_file_names ( configs ):
    """Generate intermedia file names.

    """
    sampleid = configs["sample.sample_id"]
    # decide the number of replicates as the number of comma in TREATMENT path:
    rdt = configs["data.raw_data_type"]
    # check if the raw data format matches raw data type
    if rdt == "seq":
        configs["data.number_replicates"] = len(configs["data.treatment_seq_file_path"].split(","))
    elif rdt == "alignment":
        configs["data.number_replicates"] = len(configs["data.treatment_ali_file_path"].split(","))
        
    nr = configs["data.number_replicates"]
    info("Replicates in your sample is: %d" % nr)

    # for replicates, we don't perform QC analysis, so we only keep
    # peaks.bed and treat.wig files for replicates
    configs["bowtie.treat_output_replicates"] = [sampleid+"_rep"+str(i)+"_treat.sam" for i in range(1,nr+1)]
    configs["bowtie.treat_output"] = sampleid+"_treat.sam"
    configs["bowtie.control_output"] = sampleid+"_control.sam"    

    configs["samtools.treat_output_replicates"] = [sampleid+"_rep"+str(i)+"_treat.bam" for i in range(1,nr+1)]
    configs["samtools.treat_output"] = sampleid+"_treat.bam"
    configs["samtools.control_output"] = sampleid+"_control.bam"    

    configs["macs.output_xls"] = sampleid+"_peaks.xls"
    configs["macs.output_bed_replicates"] = [sampleid+"_rep"+str(i)+"_peaks.bed" for i in range(1,nr+1)]    
    configs["macs.output_bed"] = sampleid+"_peaks.bed"
    configs["macs.output_summits"] = sampleid+"_summits.bed"
    configs["macs.output_treat_wig_replicates"] = [sampleid+"_rep"+str(i)+"_treat.wig" for i in range(1,nr+1)]        
    configs["macs.output_treat_wig"] = sampleid+"_treat.wig"
    configs["macs.output_control_wig"] = sampleid+"_control.wig"    

    configs["ceas.output_xls"] =  sampleid+"_ceas.xls"
    configs["ceas.output_pdf"] = sampleid+"_ceas.pdf"
    configs["ceas.output_R"] = sampleid+"_ceas.R"    

    configs["venn.replicates_output_png"] = sampleid+"_venn_replicates.png"
    configs["venn.dhs_output_png"] = sampleid+"_venn_dhs.png"    

    configs["correlation.output_pdf"] = sampleid+"_cor.pdf"
    configs["correlation.output_R"] = sampleid+"_cor.R"

    configs["conservation.output_bmp"] = sampleid+"_conserv.bmp"
    configs["conservation.output_R"] = sampleid+"_conserv.R"    

    configs["seqpos.output_zip"] = sampleid+"_seqpos.zip"
    return configs

# wrapper to run command 
def run_cmd ( command ):
    info ("Run: %s" % command)
    subpcall (command,shell=True)
    return

# All the pipeline command calls
def step1_bowtie (configs):
    """Step1: run bowtie to map reads to genome.
    
    """
    # check the startstep whether to pass this step

    if configs["others.startstep"] <= 1 and 1 <= configs["others.endstep"]:
        info("Step 1: BOWTIE...")
    else:
        info("Step 1 Bowtie is skipped as requested.")
        return False
    if configs["data.raw_data_type"] == "seq":
        pass
    else:
        info("Since raw data typs is not sequencing, bowtie is skipped.")
        return False

    # check the input
    tfiles = configs["data.treatment_seq_file_path"].split(",")
    cfiles = configs["data.control_seq_file_path"].split(",")
    
    for tfile in tfiles:
        if not os.path.isfile(tfile):
            error("Input for bowtie is missing: %s can't be found!" % tfile)
            sys.exit(1)
    for cfile in cfiles:
        if not os.path.isfile(cfile):
            error("Input for bowtie is missing: %s can't be found!" % cfile)
            sys.exit(1)

    idata_format = configs["data.raw_data_format"]

    # combine input data
    if len(cfiles)>=1:
        combined_input_file = "combined_input." + idata_format
        command_line = "cat "+ " ".join(cfiles) + " > " + combined_input_file
        run_cmd(command_line)
    
    # run bowtie
    if idata_format == "fastq":
        bowtie_format_option = " -q "
    elif idata_format == "fasta":
        bowtie_format_option = " -f "
    elif idata_format == "sequences":
        bowtie_format_option = " -r "

    if configs["bowtie.bowtie_max_alignment"]:
        bowtie_max_alignment_option = " -m "+configs["bowtie.bowtie_max_alignment"]+" "
    else:
        bowtie_max_alignment_option = " -m 1 "        
    # for each replicates of treatment:
    for i in range(1,configs["data.number_replicates"]+1):
        command_line = configs["bowtie.bowtie_main"]+" -S "+bowtie_format_option+bowtie_max_alignment_option+configs["bowtie.bowtie_genome_index_path"]+" "+tfiles[i-1]+" "+configs["bowtie.treat_output_replicates"][i-1]
        run_cmd(command_line)
        # convert sam to bam
        command_line = configs["samtools.samtools_main"]+" view -bt "+configs["samtools.samtools_chrom_len_path"]+" "+configs["bowtie.treat_output_replicates"][i-1]+" > "+configs["samtools.treat_output_replicates"][i-1]
        run_cmd(command_line)

    # combine replicates:
    command_line = "cat "+ " ".join(configs["bowtie.treat_output_replicates"]) + " > " + configs["bowtie.treat_output"]
    run_cmd(command_line)
    # convert sam to bam
    command_line = configs["samtools.samtools_main"]+" view -bt "+configs["samtools.samtools_chrom_len_path"]+" "+configs["bowtie.treat_output"]+" > "+configs["samtools.treat_output"]
    run_cmd(command_line)

    # for the control data:
    if len(cfiles)>=1:
        combined_input_file = "combined_input." + idata_format
        command_line = configs["bowtie.bowtie_main"]+" -S "+bowtie_format_option+bowtie_max_alignment_option+configs["bowtie.bowtie_genome_index_path"]+" "+combined_input_file+" "+configs["bowtie.control_output"]
        run_cmd(command_line)
        # convert sam to bam
        command_line = configs["samtools.samtools_main"]+" view -bt "+configs["samtools.samtools_chrom_len_path"]+" "+configs["bowtie.control_output"]+" > "+configs["samtools.control_output"]
        run_cmd(command_line)
        run_cmd("rm -f %s" % combined_input_file)
    
    return True

def step2_macs (configs):
    """Step2: run MACS to call peaks on replicates, then run again after combine replicates.
    
    """
    # check the startstep whether to pass this step

    if configs["others.startstep"] <= 2 and 2 <= configs["others.endstep"]:
        info("Step 2: MACS...")
    else:
        info("Step 2 MACS is skipped as requested.")
        return False
    if configs["data.raw_data_type"] == "alignment":
        _step2_macs_alignment(configs)
    elif configs["data.raw_data_type"] == "seq":
        _step2_macs_seq(configs)
    else:
        info("Since raw data typs is not alignment or sequences, step 2 MACS is skipped.")
        return False

def _step2_macs_alignment (configs):
    """Step2 MACS if the raw data type is alignment.
    
    """
    # check the input
    tfiles = configs["data.treatment_ali_file_path"].split(",")
    cfiles = configs["data.control_ali_file_path"].split(",")

    for tfile in tfiles:
        if not os.path.isfile(tfile):
            error("Input for MACS is missing: %s can't be found!" % tfile)
            sys.exit(1)
    if configs["data.has_control"]:
        for cfile in cfiles:
            if not os.path.isfile(cfile):
                error("Input for MACS is missing: %s can't be found!" % cfile)
                sys.exit(1)

    idata_format = configs["data.raw_data_format"]

    # combine treatment/input data
    if len(tfiles)>=1:
        if idata_format == "BAM":
            # samtools merge command
            combined_treat_ali_file = "combined_treat.bam"
            command_line = "samtools merge "+combined_treat_ali_file+" "+" ".join(tfiles)
            run_cmd(command_line)
        elif idata_format == "SAM" or idata_format == "BED":
            # samtools merge command
            combined_treat_ali_file = "combined_treat."+idata_format
            command_line = "cat "+ " ".join(tfiles) + " > " + combined_treat_ali_file
            run_cmd(command_line)
        else:
            error("Format %s not recognized!" % idata_format)
            sys.exit(1)

    if configs["data.has_control"]:
        if idata_format == "BAM":
            # samtools merge command
            combined_input_ali_file = "combined_input.bam"
            command_line = "samtools merge "+combined_input_ali_file+" "+" ".join(cfiles)
            run_cmd(command_line)
        elif idata_format == "SAM" or idata_format == "BED":
            # samtools merge command
            combined_input_ali_file = "combined_input."+idata_format
            command_line = "cat "+ " ".join(cfiles) + " > " + combined_input_ali_file
            run_cmd(command_line)
        else:
            error("Format %s not recognized!" % idata_format)
            sys.exit(1)

    # run MACS, first for each replicate
    for i in range(1,configs["data.number_replicates"]+1):
        if configs["data.has_control"]:
            # run MACS w/ control
            command_line = configs["macs.macs_main"]+" -w -S -t "+tfiles[i-1]+" -c "+ combined_input_ali_file + " -n "+configs["sample.sample_id"]+"_rep"+str(i)
            run_cmd(command_line)
            # copy out and rename the wiggle file
            command_line = "zcat "+configs["sample.sample_id"]+"_rep"+str(i)+"_MACS_wiggle/treat/"+configs["sample.sample_id"]+"_rep"+str(i)+"_treat_afterfiting_all.wig.gz > "+configs["macs.output_treat_wig_replicates"][i-1]
            run_cmd(command_line)
        else:
            # run MACS w/o control
            command_line = configs["macs.macs_main"]+" -w -S -t "+tfiles[i-1]+" -n "+configs["sample.sample_id"]+"_rep"+str(i)
            run_cmd(command_line)            
            # copy out and rename the wiggle file
            command_line = "zcat "+configs["sample.sample_id"]+"_rep"+str(i)+"_MACS_wiggle/treat/"+configs["sample.sample_id"]+"_rep"+str(i)+"_treat_afterfiting_all.wig.gz > "+configs["macs.output_treat_wig_replicates"][i-1]
            run_cmd(command_line)

    # run MACS for the combined treatment
    if configs["data.number_replicates"] == 1:
        # no need to run MACS again, simply copy the previous results
        command_line = "cp "+configs["sample.sample_id"]+"_rep1_peaks.xls"+" "+configs["macs.output_xls"]
        run_cmd(command_line)
        command_line = "cp "+configs["sample.sample_id"]+"_rep1_peaks.bed"+" "+configs["macs.output_bed"]
        run_cmd(command_line)
        command_line = "cp "+configs["sample.sample_id"]+"_rep1_summits.bed"+" "+configs["macs.output_summits"]
        run_cmd(command_line)
        command_line = "cp "+configs["macs.output_treat_wig_replicates"][0]+" "+configs["macs.output_treat_wig"]
        run_cmd(command_line)
        if configs["data.has_control"]:
            # 'copy' the wiggle file for control
            command_line = "zcat "+configs["sample.sample_id"]+"_rep"+str(i)+"_MACS_wiggle/control/"+configs["sample.sample_id"]+"_rep1_control_afterfiting_all.wig.gz > "+configs["macs.output_control_wig"]
            run_cmd(command_line)
    else:
        # run MACS on combined alignment files
        if configs["data.has_control"]:
            command_line = configs["macs.macs_main"]+" -w -S -t "+combined_treat_ali_file+" -c "+combined_input_ali_file+" -n "+configs["sample.sample_id"]
            run_cmd(command_line)
            # copy out and rename the wiggle file
            command_line = "zcat "+configs["sample.sample_id"]+"_MACS_wiggle/treat/"+configs["sample.sample_id"]+"_treat_afterfiting_all.wig.gz > "+configs["macs.output_treat_wig"]
            run_cmd(command_line)
            command_line = "zcat "+configs["sample.sample_id"]+"_MACS_wiggle/control/"+configs["sample.sample_id"]+"_control_afterfiting_all.wig.gz > "+configs["macs.output_control_wig"]
            run_cmd(command_line)
        else:
            command_line = configs["macs.macs_main"]+" -w -S -t "+combined_treat_ali_file+" -n "+configs["sample.sample_id"]
            run_cmd(command_line)
            # copy out and rename the wiggle file
            command_line = "zcat "+configs["sample.sample_id"]+"_MACS_wiggle/treat/"+configs["sample.sample_id"]+"_treat_afterfiting_all.wig.gz > "+configs["macs.output_treat_wig"]
            run_cmd(command_line)

    # rename input files by faking samtools output
    for i in xrange(len(tfiles)):
        command_line = "mv "+tfiles[i]+" "+configs["samtools.treat_output_replicates"][i]
        run_cmd(command_line)
    command_line = "mv "+combined_treat_ali_file+" "+configs["samtools.treat_output"]
    run_cmd(command_line)
    if configs["data.has_control"]:
        command_line = "mv "+combined_input_ali_file+" "+configs["samtools.control_output"]
        run_cmd(command_line)
    return True

def _step2_macs_seq (configs):
    """Step2 MACS if the raw data type is seq. So it will use the output from step1.
    
    """
    # check the input
    t_rep_files = configs["samtools.treat_output_replicates"]
    t_comb_file = configs["samtools.treat_output"]
    c_comb_file = configs["samtools.control_output"]
    macs_genome_option = " -g "+ configs["sample.species"]+" "
    
    # run MACS, first for each replicate
    for i in range(1,configs["data.number_replicates"]+1):
        if configs["data.has_control"]:
            # run MACS w/ control
            command_line = configs["macs.macs_main"]+macs_genome_option+" -w -S -t "+t_rep_files[i-1]+" -c "+ c_comb_file + " -n "+configs["sample.sample_id"]+"_rep"+str(i)
            run_cmd(command_line)
            # copy out and rename the wiggle file
            command_line = "zcat "+configs["sample.sample_id"]+"_rep"+str(i)+"_MACS_wiggle/treat/"+configs["sample.sample_id"]+"_rep"+str(i)+"_treat_afterfiting_all.wig.gz > "+configs["macs.output_treat_wig_replicates"][i-1]
            run_cmd(command_line)
        else:
            # run MACS w/o control
            command_line = configs["macs.macs_main"]+macs_genome_option+" -w -S -t "+t_rep_files[i-1]+" -n "+configs["sample.sample_id"]+"_rep"+str(i)
            run_cmd(command_line)            
            # copy out and rename the wiggle file
            command_line = "zcat "+configs["sample.sample_id"]+"_rep"+str(i)+"_MACS_wiggle/treat/"+configs["sample.sample_id"]+"_rep"+str(i)+"_treat_afterfiting_all.wig.gz > "+configs["macs.output_treat_wig_replicates"][i-1]
            run_cmd(command_line)

    # run MACS for the combined treatment
    if configs["data.number_replicates"] == 1:
        # no need to run MACS again, simply copy the previous results
        command_line = "cp "+configs["sample.sample_id"]+"_rep1_peaks.xls"+" "+configs["macs.output_xls"]
        run_cmd(command_line)
        command_line = "cp "+configs["sample.sample_id"]+"_rep1_peaks.bed"+" "+configs["macs.output_bed"]
        run_cmd(command_line)
        command_line = "cp "+configs["sample.sample_id"]+"_rep1_summits.bed"+" "+configs["macs.output_summits"]
        run_cmd(command_line)
        command_line = "cp "+configs["macs.output_treat_wig_replicates"][0]+" "+configs["macs.output_treat_wig"]
        run_cmd(command_line)
        if configs["data.has_control"]:
            command_line = "zcat "+configs["sample.sample_id"]+"_rep"+str(i)+"_MACS_wiggle/control/"+configs["sample.sample_id"]+"_rep1_control_afterfiting_all.wig.gz > "+configs["macs.output_control_wig"]
            run_cmd(command_line)
    else:
        # run MACS on combined alignment files
        if configs["data.has_control"]:
            command_line = configs["macs.macs_main"]+macs_genome_option+" -w -S -t "+t_comb_file+" -c "+c_comb_file+" -n "+configs["sample.sample_id"]
            run_cmd(command_line)
            # copy out and rename the wiggle file
            command_line = "zcat "+configs["sample.sample_id"]+"_MACS_wiggle/treat/"+configs["sample.sample_id"]+"_treat_afterfiting_all.wig.gz > "+configs["macs.output_treat_wig"]
            run_cmd(command_line)
            command_line = "zcat "+configs["sample.sample_id"]+"_MACS_wiggle/control/"+configs["sample.sample_id"]+"_control_afterfiting_all.wig.gz > "+configs["macs.output_control_wig"]
            run_cmd(command_line)
        else:
            command_line = configs["macs.macs_main"]+macs_genome_option+" -w -S -t "+t_comb_file+" -n "+configs["sample.sample_id"]
            run_cmd(command_line)
            # copy out and rename the wiggle file
            command_line = "zcat "+configs["sample.sample_id"]+"_MACS_wiggle/treat/"+configs["sample.sample_id"]+"_treat_afterfiting_all.wig.gz > "+configs["macs.output_treat_wig"]
            run_cmd(command_line)
            
    return True

def step3_ceas (configs):
    """Step3: run CEAS for the combined result.

    """
    # check the startstep whether to pass this step

    if configs["others.startstep"] <= 3 and 3 <= configs["others.endstep"]:
        info("Step 3: CEAS...")
    else:
        info("Step 3 CEAS is skipped as requested.")
        return False

    # check the input
    peak_bed_file = configs["macs.output_bed"]
    wiggle_file = configs["macs.output_treat_wig"]

    if not os.path.isfile(peak_bed_file):
        error("Input for CEAS is missing: peak file %s can't be found!" % peak_bed_file)
        sys.exit(1)
    if not os.path.isfile(wiggle_file):
        error("Input for CEAS is missing: wiggle file %s can't be found!" % wiggle_file)
        sys.exit(1)

    # run CEAS
    if configs["ceas.ceas_genetable_path"]:
        ceas_gt_option = " -g "+configs["ceas.ceas_genetable_path"]+" "
    else:
        ceas_gt_option = ""
    if configs["ceas.ceas_promoter_sizes"]:
        ceas_sizes_option = " --sizes "+configs["ceas.ceas_promoter_sizes"]+" "
    else:
        ceas_sizes_option = ""
    if configs["ceas.ceas_bipromoter_sizes"]:
        ceas_bisizes_option = " --bisizes "+configs["ceas.ceas_bipromoter_sizes"]+" "
    else:
        ceas_bisizes_option = ""
    if configs["ceas.ceas_rel_dist"]:
        ceas_rel_dist_option = " --rel-dist "+configs["ceas.ceas_rel_dist"]+" "        
    else:
        ceas_rel_dist_option = ""
    ceas_name_option = " --name "+configs["sample.sample_id"]+"_ceas "
    
    command_line = configs["ceas.ceas_main"]+ceas_name_option+ceas_gt_option+ceas_sizes_option+ceas_bisizes_option+ceas_rel_dist_option+" -b "+peak_bed_file+" -w "+wiggle_file
    run_cmd(command_line)
    
    return True

def step4_venn (configs):
    """Step4: run Venn diagram: 1) compare replicates; 2) overlap with human/mouse DHS.

    """
    # check the startstep whether to pass this step

    if configs["others.startstep"] <= 4 and 4 <= configs["others.endstep"]:
        info("Step 4: Venn diagram...")
    else:
        info("Step 4 Venn diagram is skipped as requested.")
        return False

    # check the input
    peak_bed_file_replicates = configs["macs.output_bed_replicates"]
    peak_bed_file = configs["macs.output_bed"]
    DHS_bed_file = configs["venn.dhs_bed_path"]

    for f in peak_bed_file_replicates:
        if not os.path.isfile(f):
            error("Input for Venn diagram is missing: peak file %s can't be found!" % f)
            sys.exit(1)
    if not os.path.isfile(peak_bed_file):
        error("Input for Venn diagram is missing: peak file %s can't be found!" % peak_bed_file)
        sys.exit(1)
    if not os.path.isfile(DHS_bed_file):
        error("Input for Venn diagram is missing: peak file %s can't be found!" % DHS_bed_file)
        sys.exit(1)
        
    if configs["data.number_replicates"] > 3:
        # can't process > 3 replicates. This may cause error in the following steps, so simply skip it.
        warn("Venn diagram can't process > 3 files!")
        command_line = "touch "+configs["venn.replicates_output_png"]
        run_cmd(command_line)
    elif configs["data.number_replicates"] == 1:
        # no replicates, no Venn diagram
        info("No replicates, so no Venn diagram for replicates")
        command_line = "touch "+configs["venn.replicates_output_png"]
        run_cmd(command_line)
    else:
        # run Venn diagram for replicates
        command_line = configs["venn.venn_diagram_main"]+" -t Overlap_of_Replicates"+" "+" ".join(peak_bed_file_replicates)+" "+" ".join(map(lambda x:"-l replicate_"+str(x),xrange(1,configs["data.number_replicates"]+1)))
        run_cmd(command_line)
        command_line = "mv venn_diagram.png "+configs["venn.replicates_output_png"]
        run_cmd(command_line)

    # run Venn for DHS overlap
    command_line = configs["venn.venn_diagram_main"]+" -t Overlap_with_DHS"+" -l Peaks "+peak_bed_file+" -l DHS "+DHS_bed_file
    run_cmd(command_line)
    command_line = "mv venn_diagram.png "+configs["venn.dhs_output_png"]
    run_cmd(command_line)
    
    return True

def step5_cor (configs):
    """Step5: run correlation plot tool on replicates.

    """
    # check the startstep whether to pass this step

    if configs["others.startstep"] <= 5 and 5 <= configs["others.endstep"]:
        info("Step 5: Correlation plot...")
    else:
        info("Step 5 Correlation plot is skipped as requested.")
        return False

    # check the input
    wiggle_file_replicates = configs["macs.output_treat_wig_replicates"]

    for f in wiggle_file_replicates:
        if not os.path.isfile(f):
            error("Input for correlation plot is missing: wiggle file %s can't be found!" % f)
            sys.exit(1)

    if configs["data.number_replicates"] == 1:
        # no replicates, no correlation plot
        info("No replicates, so no correlation plot for replicates")
        command_line = "touch "+configs["correlation.output_pdf"]
        run_cmd(command_line)
        command_line = "touch "+configs["correlation.output_R"]
        run_cmd(command_line)
        return True
        
    # parameters
    if configs["correlation.wig_correlation_step"]:
        cor_step_option = " -s "+configs["correlation.wig_correlation_step"]+" "
    else:
        cor_step_option = " -s 10 "

    if configs["correlation.wig_correlation_method"]:
        cor_method_option = " -m "+configs["correlation.wig_correlation_method"]+" "
    else:
        cor_method_option = " -m mean"

    if configs["correlation.wig_correlation_min"]:
        cor_min_option = " --min-score "+configs["correlation.wig_correlation_min"]+" "
    else:
        cor_min_option = " --min-score 2 "

    if configs["correlation.wig_correlation_max"]:
        cor_max_option = " --max-score "+configs["correlation.wig_correlation_max"]+" "
    else:
        cor_max_option = " --max-score 50 "

    cor_db_option = " -d "+configs["sample.assembly_name"]+" "

    # run correlation plot for replicates
    command_line = configs["correlation.wig_correlation_main"]+cor_db_option+cor_step_option+\
                   cor_method_option+cor_min_option+cor_max_option+\
                   " -r "+configs["correlation.output_R"]+" "+" ".join(wiggle_file_replicates)+\
                   " "+" ".join(map(lambda x:"-l replicate_"+str(x),xrange(1,configs["data.number_replicates"]+1)))
    run_cmd(command_line)
    command_line = "Rscript "+configs["correlation.output_R"]
    run_cmd(command_line)

    return True

def step6_conserv (configs):
    """Step6: conservation plot on peak summits of combined runs.

    """
    if configs["others.startstep"] <= 6 and 6 <= configs["others.endstep"]:
        info("Step 6: Conservation plot...")
    else:
        info("Step 6 Conservation plot is skipped as requested.")
        return False

    # check the input
    peak_summits_file = configs["macs.output_summits"]

    if not os.path.isfile(peak_summits_file):
        error("Input for conservation plot is missing: peak file %s can't be found!" % peak_summits_file)
        sys.exit(1)
        
    # run Venn diagram for replicates
    command_line = configs["conservation.conserv_plot_main"]+" -t Conservation_at_summits"+" -d "+configs["conservation.conserv_plot_phast_path"]+" -l Peak_summits "+peak_summits_file
    run_cmd(command_line)
    command_line = "mv tmp.R "+configs["conservation.output_R"]
    run_cmd(command_line)
    command_line = "mv tmp.bmp "+configs["conservation.output_bmp"]
    run_cmd(command_line)
    return True

def step7_seqpos (configs):
    """Step7: SeqPos on the top 1000 binding sites summits against TRANSFAC, JASPAR and de novo motifs.
    
    """
    if configs["others.startstep"] <= 7 and 7 <= configs["others.endstep"]:
        info("Step 7: Conservation plot...")
    else:
        info("Step 7 SeqPos is skipped as requested.")
        return False

    # check the input
    peak_summits_file = configs["macs.output_summits"]

    if not os.path.isfile(peak_summits_file):
        error("Input for SeqPos is missing: peak file %s can't be found!" % peak_summits_file)
        sys.exit(1)

    # generate top # of peaks
    psf_fhd = open(peak_summits_file)
    p_list = []
    for i in psf_fhd:
        p_list.append( (i,float(i.split()[-1])) )
    top_n = int(configs["seqpos.seqpos_top_peaks"])
    top_n_summits = map(lambda x:x[0],sorted(p_list,key=lambda x:x[1],reverse=True)[:top_n])
    top_n_summits_file = "top"+str(top_n)+"_summits.bed"
    top_n_summits_fhd = open(top_n_summits_file,"w")
    for i in top_n_summits:
        top_n_summits_fhd.write(i)
    top_n_summits_fhd.close()
    info("Top %d summits are written to %s" % (top_n,top_n_summits_file))
    
    # run SeqPos: use the current daisy version as standard, you may need to modify the options.
    # options
    if configs["seqpos.seqpos_width"]:
        seqpos_width_option = " -w "+configs["seqpos.seqpos_width"]+" "
    else:
        seqpos_width_option = " -w 600 "
    if configs["seqpos.seqpos_pvalue_cutoff"]:
        seqpos_pvalue_option = " -p "+configs["seqpos.seqpos_pvalue_cutoff"]+" "
    else:
        seqpos_pvalue_option = " -p 0.001 "
    if configs["seqpos.seqpos_motif_db_selection"]:
        seqpos_motif_option = " -m "+configs["seqpos.seqpos_motif_db_selection"]+" "
    else:
        seqpos_motif_option = " -m transfac.xml,pbm.xml,jaspar.xml "
    seqpos_filter_option = " -s "+configs["sample.species"]
    
    command_line = configs["seqpos.seqpos_main"]+" -d "+seqpos_width_option+seqpos_pvalue_option+seqpos_motif_option+seqpos_filter_option+" "+top_n_summits_file+" "+configs["sample.assembly_name"]
    run_cmd(command_line)
    command_line = "zip -r "+configs["seqpos.output_zip"]+" results/"
    run_cmd(command_line)
    return True

def step8_package_result ( configs ):
    """Package result and generate a summary file. -> a subfolder
    named as sample#sample_id and zip it as sample#sample_id.tar.bz2

    """
    if configs["others.startstep"] <= 8 and 8 <= configs["others.endstep"]:
        info("Step 8: Packaging results...")
    else:
        info("Step 8 Packaging is skipped as requested.")
        return False

    subfolder = "sample"+configs["sample.sample_id"]
    command_line = "mkdir "+subfolder
    run_cmd(command_line)
    
    all_files = []
    e = all_files.extend
    a = all_files.append
    
    e(configs["samtools.treat_output_replicates"]) # treatment bam files for replicates
    a(configs["samtools.treat_output"])            # treatment bam file for combined
    a(configs["samtools.control_output"])          # control bam file for combined
    a(configs["macs.output_xls"])
    e(configs["macs.output_bed_replicates"])
    a(configs["macs.output_bed"])
    a(configs["macs.output_summits"])
    e(configs["macs.output_treat_wig_replicates"])
    a(configs["macs.output_treat_wig"])
    a(configs["macs.output_control_wig"])
    a(configs["ceas.output_xls"])
    a(configs["ceas.output_pdf"])
    a(configs["ceas.output_R"])
    a(configs["venn.replicates_output_png"])
    a(configs["venn.dhs_output_png"])
    a(configs["correlation.output_pdf"])
    a(configs["correlation.output_R"])
    a(configs["conservation.output_bmp"])
    a(configs["conservation.output_R"])
    a(configs["seqpos.output_zip"])

    f = open(configs["macs.output_bed"])
    number_of_peaks = len(f.readlines())
    logfhd.write("----\nNumber of peaks: %d\n" % number_of_peaks)
    f.close()
    f = open(configs["macs.output_xls"])
    n = 0
    for l in f:
        l.rstrip()
        if l.startswith("# d = "):
            n = int(l[6:])
            break
    logfhd.write("Value d from MACS is %d\n----\n" % n)
    f.close()

    command_line = "mv "+" ".join(all_files)+" "+subfolder+"/"
    run_cmd(command_line)
    command_line = "cp log "+subfolder+"/"
    run_cmd(command_line)
    command_line = "tar -jcf "+subfolder+".tar.bz2 "+subfolder+"/"
    run_cmd(command_line)


# ------------------------------------
# Main function
# ------------------------------------
def main():
    usage = "usage: %prog <config file>"
    description = "ChIP-seq pipeline demo"
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-s","--step",dest="startstep",type="int",
                         help="Select the step to start from, use it if you want to resume the pipeline after any interruption. Step 1: bowtie mapping; step 2: macs peak calling; step 3: CEAS analysis; step 4: Venn diagram; step 5: correlation; step 6: conservation; step 7: seqpos; step 8: packaging the results. Default is blank, so that the pipeline will start from the step 1.", default=1)
    optparser.add_option("-S",dest="stopatstep",action="store_true",default=False,
                         help="If True, the pipeline will stop after it finished the step specified by -s option.")    
    (options,args) = optparser.parse_args()

    startstep = options.startstep
    stopatstep = options.stopatstep

    if not args or startstep<1 or startstep>8:
        optparser.print_help()
        sys.exit(1)

    # parse the config file
    config_file = args[0]
    configs = read_config(config_file)

    # decide the start step and end step of the pipeline
    configs["others.startstep"] = startstep
    if stopatstep:
        configs["others.endstep"] = startstep
    else:
        configs["others.endstep"] = 8

    configs["others.cwd"] = os.getcwd()
    
    if not check_conf_validity(configs):
        error("Exit as failing config file validity check!")
        sys.exit(1)
    else:
        info("Pass the conf file validity check!")

    if not check_file_validity(configs):
        error("Exit as failing data file validity check!")
        sys.exit(1)
    else:
        info("Pass the data file validity check!")

    if not check_cmd_validity(configs):
        error("Exit as failing command validity check!")
        sys.exit(1)
    else:
        info("Pass the command line validity check!")

    if not check_lib_validity(configs):
        error("Exit as failing tool libraries validity check!")
        sys.exit(1)
    else:
        info("Pass the tool libraries file validity check!")

    info("All checks are passed! Prepare the pipeline...")

    steps = [step1_bowtie,step2_macs,step3_ceas,step4_venn,step5_cor,step6_conserv,step7_seqpos,step8_package_result]
    prepare_output_file_names(configs)
    steps[0](configs)
    steps[1](configs)
    steps[2](configs)
    steps[3](configs)
    steps[4](configs)
    steps[5](configs)                    
    steps[6](configs)
    steps[7](configs)    
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)

    
logfhd.close()
