#!/usr/bin/env python
# Time-stamp: <2011-05-05 18:31:45 Tao Liu>

"""Module Description

Copyright (c) 2011 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: 0.1
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------

import sys
from optparse import OptionParser
 
from Bio import Entrez
 
# ------------------------------------
# constants
# ------------------------------------
# *Always* tell NCBI who you are
Entrez.email = "your email here"

# ------------------------------------
# Misc functions
# ------------------------------------
def search_genes(id_list,search_field):
    """Use ESearch to convert RefSeq or Gene symbols to standard
    Entrez IDs.

    A request to esearch.cgi is like:
    http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=ID_LIST[SEARCH_FIELD]

    Return a list of Entrez IDs.
    """
    term = " OR ".join(map(lambda x:x+"["+search_field+"]",id_list))
    esearch_result = Entrez.esearch(db="gene",term=term,retmod="xml")
    parsed_result = Entrez.read(esearch_result)
    return parsed_result['IdList']

def fetch_genes(id_list):
    """Fetch Entrez Gene records using Bio.Entrez, in particular epost
    (to submit the data to NCBI) and efetch to retrieve the
    information, then use Entrez.read to parse the data.

    Returns a list of parsed gene records.
    """
 
    request = Entrez.epost("gene",id=",".join(id_list))
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        print "An error occurred while retrieving the annotations."
        print "The error returned was %s" % e
        sys.exit(-1)
 
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    efetch_result = Entrez.efetch(db="gene", webenv=webEnv, query_key = queryKey, retmode="xml")
    genes = Entrez.read(efetch_result)
    #print "Retrieved %d records for %d genes" % (len(genes),len(id_list))
    return genes

def parse_genes(genes):
    """Parse various gene information including:

    1. Species name (taxonomy name)
    2. Entrez gene ID
    3. Official symbol
    4. RefSeq IDs
    5. Offical full name

    Basically, just to go through the parsed xml data.... A big headache to figure it out...

    Return a list of dictionary.
    """
    gene_info_list = []
    for gene_data in genes:
        gene_info = {}
        # get entrez ID
        try:
            gene_info["entrez_id"] = gene_data["Entrezgene_track-info"]["Gene-track"]["Gene-track_geneid"]
        except KeyError:
            gene_info["entrez_id"] = ""
            continue
        gene_info["refseq_ids"] = []
        for comment in gene_data.get("Entrezgene_comments",[]):
            # look for refSeq annotation
            if comment.get("Gene-commentary_heading",None) == "NCBI Reference Sequences (RefSeq)":
                # get sub-comments
                for subcomment in comment.get("Gene-commentary_comment",[]):
                    for product in subcomment.get("Gene-commentary_products",[]):
                        if product.get("Gene-commentary_heading",None) == "mRNA Sequence":
                            gene_info["refseq_ids"].append(product.get("Gene-commentary_accession",""))
        # get properties
        gene_info["official_symbol"] = "" # optional
        gene_info["official_full_name"] = "" # optional
        for gene_property in gene_data.get("Entrezgene_properties",[]):
            if gene_property.get("Gene-commentary_label",None) == "Nomenclature":
                for sub_property in gene_property["Gene-commentary_properties"]:
                    if sub_property.get("Gene-commentary_label",None)  == "Official Symbol":
                        gene_info["official_symbol"] = sub_property.get("Gene-commentary_text","")
                    if sub_property.get("Gene-commentary_label",None)  == "Official Full Name":
                        gene_info["official_full_name"] = sub_property.get("Gene-commentary_text","")

        # get taxname
        try:
            gene_info["taxname"] = gene_data["Entrezgene_source"]["BioSource"]["BioSource_org"]["Org-ref"]["Org-ref_taxname"]
        except KeyError:
            gene_info["taxname"] = ""
            continue
        gene_info_list.append(gene_info)

    return gene_info_list

def print_genes (gene_info_list):
    """Print out parsed entrez gene information in tab-delimited way.

    """
    # header
    print "%s\t%s\t%s\t%s\t%s" % ("TaxonomyName","EntrezID","OfficialSymbol","RefSeqIDs","OfficialFullName")
    for g in gene_info_list:
        print "%s\t%s\t%s\t%s\t%s" % (g["taxname"],g["entrez_id"],g["official_symbol"],",".join(g["refseq_ids"]),g["official_full_name"])

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def main():
    usage = "usage: %prog [options]"
    description = "Use NCBI web API to convert gene ids between different identifier types."
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-i","--id",dest="ids",type="string",action="append",
                         help="Gene id, according to identifier setting of input, can be Entrez, RefSeq, or Gene symbol. Multiple ids are allowed.")
    optparser.add_option("-a","--itype",dest="itype",default="entrez",
                         help="Identifier type of your input ids. Can be 'entrez', 'refseq', or 'symbol'. Default: 'entrez'.")
    (options,args) = optparser.parse_args()
    if not options.ids:
        optparser.print_help()
        sys.exit(-1)
    input_id_list = options.ids
    if options.itype == "refseq":
        entrez_id_list = search_genes(input_id_list,"ACCN")
    elif options.itype == "symbol":
        entrez_id_list = search_genes(input_id_list,"GENE")
    elif options.itype == "entrez":
        entrez_id_list = input_id_list
    
    entrez_id_genes = fetch_genes(entrez_id_list)
    parsed_genes = parse_genes(entrez_id_genes)
    print_genes(parsed_genes)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
