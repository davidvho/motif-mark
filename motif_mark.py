#!/usr/bin/env python

##################################
## ARGPARSE ######################
##################################

import argparse

def mm():
	parser = argparse.ArgumentParser(description='A script where given a FASTA file where the sequence of each record is the form intronEXONintron, script will output a visual of where motifs are located in each sequence. The output is a PNG and SVG in the working directory. Requires the following python packages: re, SeqIO, cairo, random, numpy, colosys. Last updated: March 18, 2018.')
	parser.add_argument("-f", "--fasta", help="Required. FASTA file where the sequence of each record is the form intronEXONintron.", required=True, type=str)
	parser.add_argument("-m", "--motifs", help="Required. One motif per line. 10 motifs max currently.", required=True, type=str)
	return parser.parse_args()

args = mm()

# File
fasta = args.fasta
motif_file = args.motifs

########################################
## DEFINITIONS AND FUNCTIONS ###########
########################################

## 10: red, blue, green, orange, magenta, chocolate, purple, light steel blue, lawn green, tomato
colors = (255,0,0), (0,0,255), 	(124,252,0), (255,165,0),(255,0,255), (138,43,226), (210,105,30), (176,196,222),(124,252,0) ,  	(255,99,71)


iupac = {"A": "[Aa]", "T":"[TtUu]", "C":"[Cc]", "G":"[Gg]", "U":"[UuTt]",
              "R": "[AaGg]", "Y":"[TtCcUu]", "S":"CcGg", "W":"AaTtUu",
              "K":"[GgTtUu]", "M":"[AaCc]", "B":"[CcGgTtUu]", "D":"[AaGgTtUu]",
              "H":"[AaCcTtUu]", "V":"[AaCcGg]", "N":"[AaTtCcGgUu]"}


def motifs_regex(file):
    '''Takes in a text file of motifs, one for every line and convert them to regex-able strings'''
    motifs_regex_list = []
    with open(file, "r") as fh:
        for line in fh:
            motif = ""
            line = line.strip("\n").upper()
            for char in line:
                motif += iupac[char]
            motifs_regex_list.append(motif)
                
    return motifs_regex_list


def motifs(file):
    '''Takes in a text file of motifs, and makes a dictionary with regex : motif'''
    motifs = {}
    motifs_regexx = []
    with open(file, "r") as fh:
        for line in fh:
            motif = ""
            line = line.strip("\n").upper()
            for char in line:
                motif += iupac[char]
            motifs_regexx.append(motif)
            motifs[motif] = line
                
    return motifs


def len_gene(file):
    '''Given a FASTA file, outputs the sequence length for each record in a dictionary'''
    length = {}
    with open(fasta) as fp_in:
        for record in SeqIO.parse(fp_in, "fasta"):
            record = ("{},{}".format(record.id, record.seq))
            gene_name = record.split(",")[0]
            seq = record.split(",")[1]
            length[gene_name] = len(seq)
    return length   


def draw_gene_exon(length, exon_start, exon_length, y):
    '''given a length of the gene and the exon start and exon length, draw it'''

    start = 100
    exon = start + exon_start
    exon2 = start + exon_length
    
    context.set_line_width(3)
    context.move_to(start,y)        #(x,y)
    context.line_to(length+start,y)
    context.stroke()
    
    context.set_line_width(12)
    context.move_to(exon,y)        #(x,y)
    context.line_to(exon2,y)
    context.stroke()
    

##########################
## MAIN ##################
##########################

import re
from Bio import SeqIO
import cairo
import random
import numpy as np
import colorsys

## gene lengths
gene_length = len_gene(motifs)

## go through the fasta file and get dictionaries
record_d = {}  ## gene_name: [intron, exon, intron]
exon_loc = {}  ## gene_name : exon_start pos, length
with open(fasta) as fp_in:
    for record in SeqIO.parse(fp_in, "fasta"):


        record = ("{},{}".format(record.id, record.seq))
        
        gene_name = record.split(",")[0]
        seq = record.split(",")[1]
        
        seq_split = []
        exon = re.search('[A-Z]', seq)
        
        seq_split.append(seq[0:exon.start()])
        
        intron2 = re.search('[a-z]', seq[exon.start():])
        
        seq_split.append(seq[exon.start():exon.start()+intron2.start()])
        seq_split.append(seq[exon.start()+intron2.start():])

        record_d[gene_name] = seq_split
        exon_loc[gene_name] = exon.start(), exon.start()+intron2.start()


### go through the motifs regex and save the following information in dictionary--  
### GENE_MOTIF : position.found, position.found, ...
where = {}
motifs_regex_list = motifs_regex(motif_file)
og_motifs = motifs(motif_file)

for gene in record_d.keys():
    seq = "".join(record_d[gene])
    for motif in motifs_regex_list:
        loc = re.finditer(motif, seq)
        where[gene+"_"+og_motifs[motif]] = []
        for l in loc:
            where[gene+"_"+og_motifs[motif]].append(l.start())

## assign colors to motifs
c_d = {}
with open(motif_file, "r") as fh:
        no = 0
        for line in fh:
            motif = ""
            line = line.strip("\n").upper()
            c_d[line] = colors[no]
            no += 1
            
## max gene length to scale everything
longest = gene_length[max(gene_length, key=lambda i: gene_length[i])]


## start drawing!
surface = cairo.SVGSurface("output.svg", longest*1.2, len(record_d)*60+50)   ## width x height, use longest gene as setting for width, use no of genes as setting for height
context = cairo.Context(surface)

## for each gene, draw the gene proportional to length, and the exon, label the gene to the left
p = 25
right = 100
for i in gene_length.keys(): 
    context.set_source_rgb(0,0,0)
    context.move_to(50, p+4)
    context.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    context.set_font_size(15)
    context.show_text(i)
    draw_gene_exon(gene_length[i], exon_loc[i][0], exon_loc[i][1], p)
    
    ## put the motifs on the gene, color coded
    for og in og_motifs.values():
        m = i+"_"+og
        context.set_source_rgb(c_d[og][0]/265, c_d[og][1]/265, c_d[og][2]/265)
        for a in range(len(where[m])):
            context.set_line_width(12)
            context.move_to(where[m][a]+100,p)
            context.line_to(where[m][a]+103,p)
            context.stroke()
            
    ## each gene is 65 pixels offset lower
    p += 65

## put the legend at the bottom    
for og in og_motifs.values():
        m = i+"_"+og
        
        context.set_source_rgb(c_d[og][0]/265, c_d[og][1]/265, c_d[og][2]/265)
        
        context.move_to(right, len(gene_length)*65+20)
        context.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        #context.set_font_size(12)
        context.show_text(og)
        right += 100    

## write out a png and svg
surface.write_to_png("output.png")
surface.finish()

print("Done. See output files.")