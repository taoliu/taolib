

"""Module Description

Copyright (c) 2008 H. Gene Shin <shin@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  H. Gene Shin
@contact: shin@jimmy.harvard.edu

This module contains functions that produce R scripts corresponding to certain R functions. 
For example, R.plot(x,y) produces a R script of plot(x,y).

For the details on the parameters of the below functions, refer to R documents.
"""

# ------------------------------------
# Python modules
# ------------------------------------
import sys
from array import *
import Cistrome.Assoc.graphics as graphics


# ------------------------------------
# function - statistics
# ------------------------------------ 
def pbinom(q,size,p,lower_tail=True,log_p=False):
    """produce a script of pbinom in R"""
    
    rscript='pbinom('
    if type(q).__name__=='list' or 'tuple':
        rscript+='c(%s),' %str(q)[1:-1]
    else:
        rscript+='%s,' %(str(q))
    if type(size).__name__=='list' or 'tuple':
        rscript+='c(%s),' %str(q)[1:-1]
    else:
        rscript+='%s,' %(str(q))
    if type(p).__name__=='list' or 'tuple':
        rscript+='c(%s),' %str(q)[1:-1]
    else:
        rscript+='%s,' %(str(q))

    if not lower_tail:
        rscript+='lower.tail=FALSE,'
    if log_p:
        rscript+='log.p=TRUE,'
    rscript.rstrip(',')
    rscript+=')\n'
    
    return rscript

# ------------------------------------
# function - general
# ------------------------------------
def vector(x, varname='x'):
    """Return an one-line R script of a vector
    
    Parameters:
    1. x: a list or an array of numbers or string
    2. type: the type of elements of x
    3. varname: the name of the vector
    """
    
    return "%s <- c(%s)\n" %(varname, str(x)[1:-1])


def matrix(mat,nrow=1,ncol=1,byrow=False):
    """Given a two dimensional array, write the array in a matrix form"""
    
    nr=len(mat)
    rscript='m<-matrix(data=c('
    try:
        nc=len(mat[0])
        for m in mat:
            rscript+=str(m)[1:-1]+ ', '   
        rscript=rscript[:-2]+'), nrow=%d, ncol=%d, byrow=TRUE,' %(nr,nc)
    except TypeError:
        rscript+=str(mat)[1:-1]+','
        rscript=rscript[:-1]+'), nrow=%d, ncol=%d,' %(nrow,ncol)
        if byrow: rscript+='byrow=TRUE,'
    
    rscript=rscript[:-1]+')\n'
    
    return rscript

# ------------------------------------
# function - graphhics
# ------------------------------------   
def barplot(height,names=None,beside=False,horiz=False,col=None,main=None,xlab=None,\
            ylab=None,xlim=None,ylim=None,border=True,space=None, cex_names=None, cex_axis=None):
    
    # convert the data into R script
    if not height: rscript='height<-c()\n'
    else:
        if type(height[0])==list:
            rscript=''
            for i in range(0,len(height)):
                rscript+='r%d<-c(%s)\n' %(i,str(height[i])[1:-1])
            rscript+='height<-rbind('
            for i in range(0,len(height)):
                rscript+='r%d,' %i
            rscript=rscript[:-1]+')\n'
        else:
            if type(height)==array:
                rscript='height<-c(%s)\n' %str(list(height))[1:-1]
            else:
                rscript='height<-c(%s)\n' %(str(height)[1:-1])
    if names:
        rscript+='names=c('
        for n in names:
            rscript+='"%s",' %str(n)
        rscript=rscript[:-1]+')\n'
        
   
    rscript+='mp<-barplot(height=height,names=names,'
    
    if beside: rscript+='beside=TRUE,'
    else: rscript+='beside=FALSE,'
    if horiz: rscript+='horiz=TRUE,'
    else: rscript+='horiz=FALSE,'
    if col:
        rscript+='col=c('
        for c in col: 
            rscript+='"%s",' %c
        rscript=rscript[:-1]+'),'
    if main: rscript+='main="%s",' %main
    if xlab: rscript+='xlab="%s",' %xlab
    if ylab: rscript+='ylab="%s",' %ylab
    if border: rscript+='border=TRUE,'
    else: rscript+='border=NA,'
    if space: rscript+='space=%f,' %space
    if xlim:    
        rscript+='xlim=c('
        for x in xlim:
            rscript+='%f,' %x
        rscript=rscript[:-1]+'),'
    if ylim:
        rscript+='ylim=c('
        for y in ylim:
            rscript+='%f,' %y
        rscript=rscript[:-1]+'),'
        
    if cex_names:
        rscript +='cex.names=%s,' %str(cex_names)
        
    if cex_axis:
        rscript +='cex.axis=%s,' %str(cex_axis)
    
    rscript=rscript[:-1]+')\n'
    return rscript

def pie(x,labels=None, main=None, col=None, clockwise=False, border=False, init_angle = None, cex=1, radius = 0.8):
    
    # data values
    rscript='x<-c('
    if x:
        for i in x:
            rscript+='%f,' %i
        rscript=rscript[:-1]+')\n'
    else:
        rscript+=')\n'
    
    # graphic parameters
    rscript+='pie(x=x,'
    if labels:
        if type(labels) == list or type(labels) == tuple:
            rscript+='labels=c('
            for l in labels:
                rscript+='"%s",' %l
            rscript=rscript[:-1]+'),'
        elif type(labels) == str:
            rscript += 'labels=%s,' %labels
            
    if main: rscript+='main="%s",' %main
    if col: 
        if type(col) ==str:
            rscript += 'col=%s,' %col
        elif type(col) == list or type(col) == tuple:
            rscript+='col=c('
            for c in col: 
                rscript+='"%s",' %c
            rscript=rscript[:-1]+'),'
    if clockwise: rscript+='clockwise=TRUE,'
    else: rscript+='clockwise=FALSE,'
    if border: rscript+='border=TRUE,'
    else: rscript+='border=FALSE,'
    
    if radius !=0.8:
        rscript += 'radius=%s,' %str(radius)
    
    if cex !=1:
        rscript += 'cex=%s,' %str(cex)
        
    if init_angle !=None:
        rscript += 'init.angle=%s,' %str(init_angle)
    
    rscript=rscript[:-1]+')\n'
    return rscript

def plot(x,y,tp="l", col=None, main=None, xlab=None, ylab=None, xlim=None, ylim=None, frame=True, axes=True, xaxt='s', yaxt='s', lwd=None, cex=None):
    
    if type(x)==list or type(x)==tuple or type(x)==array:
        if x:
            rscript='x<-c('
            for i in x:
                if i==i: rscript+='%f,' %i
                else: rscript+='NaN,'
            rscript=rscript[:-1]+')\n'
        else:
            rscript = 'x<-c()\n'
    elif type(x)==int or type(x)==float:
        if x==x: rscript='x<-c(%f)\n' %x
        else: rscript='x<-c(NaN)\n'
    elif type(x)==str:
        rscript='x <- %s\n' %x
    
    if type(y)==list or type(y)==tuple or type(x)==array:
        if y:
            rscript+='y<-c('
            for i in y:
                if i==i:rscript+='%f,' %i
                else: rscript+='NaN,'
            rscript=rscript[:-1]+')\n'
        else:
            rscript+='y<-c()\n'
    elif type(y)==int or type(y)==float:
        if y==y: rscript+='y<-c(%f)\n' %x
        else: rscript+='y<-c(NaN)\n'
    elif type(x)==str and y:
        rscript='y <- %s\n' %y
    
    if type(y)==str and not y:
        rscript+='plot(x, '
    else:
        rscript+='plot(x, y,'

    rscript+='type="%s",' %tp
    if main!=None: rscript+='main="%s",' %main
    if xlab!=None: rscript+='xlab="%s",' %xlab
    if ylab!=None: rscript+='ylab="%s",' %ylab
    if col:
        if type(col)==list or type(col)==tuple:
            rscript+='col=c('
            for c in col: 
                rscript+='"%s",' %c
            rscript=rscript[:-1]+'),'
        elif type(col)==str:
            rscript +='col=%s,' %col
    
    if xlim:
        rscript+='xlim=c('
        for xl in xlim:
            rscript+='%f,' %xl
        rscript=rscript[:-1]+'),'
        
    if ylim: 
        rscript+='ylim=c(' 
        for yl in ylim:
            rscript+='%f,' %yl
        rscript=rscript[:-1]+'),'
        
    if not frame:
        rscript+='frame=FALSE,'
        
    if not axes:
        rscript+='axes=FALSE,'
        
    if xaxt != 'r':
        rscript+='xaxt="%s",' %xaxt
    
    if yaxt != 'r':
        rscript+='yaxt="%s",' %yaxt
        
    if lwd:
        rscript+='lwd=%d,' %lwd
        
    if cex:
        rscript+='cex=%s,' %str(cex)
    
    rscript=rscript[:-1]+')\n'
    
    return rscript

def polygon(x, y = None, density=None, angle=45, border = None, col = None):
    """polygon function"""
    
    if type(x)==list or type(x)==tuple:
        rscript='x<-c('
        for i in x:
            if i==i: rscript+='%f,' %i
            else: rscript+='NaN,'
        rscript=rscript[:-1]+')\n'
    elif type(x)==int or type(x)==float:
        if x==x: rscript='x<-c(%f)\n' %x
        else: rscript='x<-c(NaN)\n'
    
    if type(y)==list or type(y)==tuple:
        rscript+='y<-c('
        for i in y:
            if i==i:rscript+='%f,' %i
            else: rscript+='NaN,'
        rscript=rscript[:-1]+')\n'
    elif type(y)==int or type(y)==float:
        if y==y: rscript='y<-c(%f)\n' %x
        else: rscript='y<-c(NaN)\n'
    
    rscript+='polygon(x,y,'
    if density:
        if type(density)==tuple or type(density)==list:
            rscript+='density=c('
            for d in density:
                rscript+='%d,' %d
            rscript+='),'
    
    if angle != 45:
        rscript+='angle=%d,' %angle
    
    if border:
        rscript+='border="%s",' %border
    
    if col:
        if type(col) ==list or type(col)==tuple:
            rscript += 'col=c('
            for c in col:
                rscript+= '"%s",' %str(c)
            rscript = rscript[:-1]+'),'
        elif type(col)==str:
            rscript += 'col=%s,' %str(col)
    
    rscript=rscript[:-1]+')\n'
    
    return rscript
        
        

def lines(x,y,tp="l",col=None, xlim=None, ylim=None, lwd=None):
    
    if type(x)==list or type(x)==tuple or type(x)==array:
        if x:
            rscript='x <- c('
            for i in x:
                if i==i: rscript+='%f,' %i
                else: rscript+='NaN,'
            rscript=rscript[:-1]+')\n'
        else:
            rscript='x<-c()\n'
    elif type(x)==int or type(x)==float:
        if x==x: rscript='x <-c (%f)\n' %x
        else: rscript='x <- c(NaN)\n' 
    elif type(x)==str:
        rscript='x <- %s\n' %x
    
    if type(y)==list or type(y)==tuple or type(y)==array:
        if y:
            rscript+='y<-c('
            for i in y:
                if i==i:rscript+='%f,' %i
                else: rscript+='NaN,'
            rscript=rscript[:-1]+')\n'
        else:
            rscript+='y<-c()\n'
    elif type(y)==int or type(y)==float:
        if y==y: rscript='y<-c(%f)\n' %x
        else: rscript='y<-c(NaN)\n'
    elif type(y)==str and y:
        rscript='y <- %s\n' %y
    
    if type(y)==str and not y:
        rscript+='lines(x, '
    else:
        rscript+='lines(x, y,'
    
    if xlim:
        rscript+='xlim=c(%s),' %str(list(xlim))[1:-1]
    
    if ylim:
        rscript+='ylim=c(%s),' %str(list(ylim))[1:-1]
        
    rscript+='type="%s",' %tp
    if col:
        if type(col)==list or type(col)==tuple:
            rscript+='col=c('
            for c in col: 
                rscript+='"%s",' %c
            rscript=rscript[:-1]+'),'
        elif type(col)==str:
            rscript+='col=%s,' %col
    if lwd:
        rscript+='lwd=%d,' %lwd
    
    rscript=rscript[:-1]+')\n'
    return rscript

def abline(a=None,b=None,h=None,v=None,lty=None,col=None):
    
    rscript='abline('
    if a!=None:
        rscript+='a=%f,' %a
    if b!=None:
        rscript+='b=%f,' %b
    if h!=None:
        rscript+='h=%f,' %h
    if v!=None:
        rscript+='v=%f,' %v
    if lty!=None:
        rscript+='lty=%d,' %lty
    if type(col).__name__=='list' or type(col).__name__=='tuple':
        rscript+='col=c('
        for c in col:
            rscript+='"%s",' %c
        rscript=rscript[:-1]+'),'
    elif type(col).__name__=='str':
        rscript='col="%s",' %col
    
    rscript=rscript[:-1]+')\n'
    
    return rscript

def hist(x, breaks, freq=False, right = True, density=None, main='', xlim=None, ylim=None, xlab='', ylab='', axes=True, col=None, border=None, returnClass=False):
    """Histogram"""
    
    if type(x)==list or type(x)==tuple or type(x)==array:
        rscript = 'x <- c(%s)\n' %str(list(x))[1:-1]
    elif type(x)==str and x:
        rscript = 'x <- %s\n' %x
    
    if type(breaks)==list or type(breaks)==tuple or type(x)==array:
        rscript += 'breaks <- c(%s)\n' %str(list(breaks))[1:-1]
    elif type(breaks)==int:
        rscript += 'breaks <- %d\n' %breaks
    elif type(breaks)==str and breaks:
        rscript += 'breaks <- "%s"\n' %breaks
        
    # draw or just return histogram class?
    if returnClass:
        rscript +='hs <- hist(x, breaks, '
    else:
        rscript +='hist(x, breaks, '
    
    if freq:
        rscript += 'freq=TRUE, probability=FALSE, '
    else:
        rscript += 'freq=FALSE, probability=TRUE, '
    
    if not right:
        rscript += 'right=FALSE, '
    
    if density:
        rscript += 'density=%f, ' %density
    
    if main:
        rscript += 'main="%s", ' %main
    
    if xlim:
        rscript += 'xlim=c(%s), ' %str(list(xlim))[1:-1]
    
    if ylim:
        rscript += 'ylim=c(%s), ' %str(list(ylim))[1:-1]
    
    if xlab:
        rscript += 'xlab="%s", ' %xlab
    
    if ylab:
        rscript += 'ylab="%s", ' %ylab
    
    if not axes:
        rscript += 'axes=FALSE, '
        
    if col:
        rscript += 'col=%s, ' %col
    
    if border:
        rscript += 'border=%s, ' %border
        
    rscript =rscript[:-2] + ')\n'
    
    return rscript

def seq(fr, to, by=1, rval = 's'):
    """Seq function"""
    
    rscript = ''
    if rval:
        if by != 1:
            rscript += '%s <- seq(from=%s, to=%s, by=%s)\n' %(rval, str(fr), str(to), str(by))
        else:
            rscript += '%s <- seq(from=%s, to=%s)\n' %(rval, str(fr), str(to))
    else:
        if by != 1:
            rscript += 'seq(from=%s, to=%s, by=%s)\n' %(str(fr), str(to), str(by))
        else:
            rscript += 'seq(from=%s, to=%s)\n' %(str(fr), str(to))
    
    return rscript


def legend(x,legend,pch,y=None,col=None,bty='o'):
    
    rscript='legend('
    if type(x).__name__=='str':
        rscript+='"%s",' %x
    elif type(x).__name__=='int' or type(x).__name__=='float':
        rscript+='%f,' %x
    if y: rscript+='%f,' %y
    rscript+='legend=c('
    for l in legend:
        rscript+='"%s",' %l
    rscript=rscript[:-1]+'),'
    if col:
        if type(col)==list or type(col)==tuple:
            rscript+='col=c('
            for c in col: 
                rscript+='"%s",' %c
            rscript=rscript[:-1]+'),'
        elif type(col)==str:
            rscript+='col=%s,' %col
    rscript+='pch=%d,' %pch
    rscript+='bty="%s",' %bty 
    rscript=rscript[:-1]+')\n'
    
    return rscript

        
def pdf(filename,height,width):
    
    rscript='pdf("%s",height=%.1f,width=%.1f)\n' %(filename,height,width)
    return rscript

def devoff():
    
    return 'dev.off()\n'

###
# par is underconstruction. It just adjusts margins and multiple plots in one figure currently.

def par(mfrow=None,mfcol=None,mar=None,oma=None,xpd=None,xaxt=None,yaxt=None):
    
    rscript='par('
    if mfrow: rscript+='mfrow=c('+str(mfrow)[1:-1]+'),'
    if mfcol: rscript+='mfcol=c('+str(mfcol)[1:-1]+'),'
    if mar: rscript+='mar=c('+str(mar)[1:-1]+'),'
    if oma: rscript+='oma=c('+str(oma)[1:-1]+'),'
    if xpd: rscript+='xpd=TRUE,'
    if xaxt: rscript+'xaxt="%s",' %xaxt
    if yaxt: rscript+='yaxt="%s",' %yaxt
    
    rscript=rscript[:-1]+')\n'
    return rscript

def mtext(text,side=3,line=0,at=None,font=None,outer=False,cex=None):
    
    if not text: return ''
    rscript='mtext("%s",' %text
    rscript+='side=%d,' %side
    rscript+='line=%d,' %line
    rscript+='outer=%s,' %(not outer and 'FALSE' or 'TRUE')
    
    if type(font).__name__=='int' or type(font).__name__=='float': rscript+='font=%d,' %font
    elif type(font).__name__=='tuple' or type(font).__name__=='list': 
        rscript+='font=c('
        for f in font:
            rscript+='%d,' %f
        rscript=rscript[:-1]+'),'
    
    if type(cex).__name__=='int' or type(cex).__name__=='float': rscript+='cex=%s,' %str(cex)
    elif type(cex).__name__=='tuple' or type(cex).__name__=='list': 
        rscript+='cex=c(%s),' %str(cex)[1:-1]
   
    if at:
        rscript += 'at=%s,' %str(at)
    
    rscript=rscript[:-1]+')\n'
     
    return rscript

def text(x,y,label,pos=None,offset=0.5,cex=None):
    
    if type(x).__name__=='list':
        rscript='text(x=c('+str(x)[1:-1]+'),'
    elif type(x).__name__=='str':
        rscript='text(x=%s,' %x
    elif type(x).__name__=='int' or type(x).__name__=='float':
        rscript='text(x=%.1f,' %x
    
    if type(y).__name__=='list':
        rscript+='y=c('+str(y)[1:-1]+'),'
    elif type(y).__name__=='str':
        rscript+='y=%s,' %y
    elif type(y).__name__=='int' or type(y).__name__=='float':
       rscript+='y=%.1f,' %y
    
    if type(label).__name__=='list':
        rscript+='label=c('
        for l in label:
            rscript+='"%s",' %l
        rscript=rscript[:-1]+'),'
    elif type(label).__name__=='str':
        rscript+='label="%s",' %label
        
    if pos:
        rscript+='pos=%d,' %pos
    
    if offset!=0.5:
        rscript+='offset=%.1f,' %offset
    
    if type(cex).__name__=='int' or type(cex).__name__=='float': 
        rscript+='cex=%s,' %str(cex)
    elif type(cex).__name__=='tuple' or type(cex).__name__=='list': 
        rscript+='cex=c(%s),' %str(cex)[1:-1]

    rscript=rscript[:-1]+')\n'    
    return rscript


def layout(mat,widths=None,heights=None):
    """layout"""
    
    ncol=len(mat[0])
    nrow=len(mat)
    arr=[]
    map(lambda m: arr.extend(m),mat)
    rscript='layout(matrix(c(%s), %d, %d, byrow = TRUE),' %(str(arr)[1:-1],nrow,ncol)
    if widths:
        rscript+='widths=c(%s),' %(str(widths)[1:-1])
    if heights:
        rscript+='heights=c(%s),' %(str(heights)[1:-1])
    
    rscript=rscript[:-1]+')\n'

    return rscript
     
        
def comment(message):
    
    rscript='# ' +message
    rscript+='\n'
    
    return rscript


def rainbow(n, s = 1, v = 1, start = 1, rval='cols'):
    """Return n rainbow colors
    """
    rscript = ''
    if not n: return rscript
    
    rscript += '%s <- rainbow(%d,' %(rval, n)
    if s != 1:
        rscript += 's=%s,' %str(s)
    if v != 1:
        rscript += 'v=%s,' %str(v)
    if start != 1:
        rscript += 'start=%s,' %str(start)
    
    rscript = rscript[:-1] +')\n'
    
    return rscript


def colorRampPalette(colors, bias=1, space='rgb', interpolate='linear', rval='colorPalette'):
    """ColorRampPalette function"""

    rscript = ''
    if colors:
        if type(colors) == list or type(colors) == tuple:
            rscript += '%s <- colorRampPalette(c(%s),' %(rval, str(colors)[1:-1])
        elif type(colors) == str:
            rscript += '%s <- colorRampPalette(%s,' %(rval, colors)
    
    if bias != 1:
        rscript += 'bias = %s,' %str(bias)
    
    if interpolate != 'linear':
        rscript += 'interpolate = "%s",' %str(interpolate)
    
    if rscript:
        rscript = rscript[:-1] +')\n'
    
    return rscript


def heatmap_bar(x, y, fromto, ylim, bot=0, top=2, vertical=False, cmap='cmap'):
    """Draw a bar of heatmap.
    
    Parameters;
    1. x: x values
    2. y: y values
    3. fromto: the first and last points of the bar. If the first point is smaller than x's smallest, between it and the first of x there is no color
    4. ylim: ylimits
    5. bot: bottom of the bar
    6. top: top of the bar
    7. vertical: the bar stands vertical or lies horizontal?
    8. cmap: the colormap variable (vector) name in the R script. WARNING: to use heatmap_bar, a colormap must be given prior to this function and correctly referred to. 
    """
 
    # quantize y values
    #qy = graphics.quantize(list(y), ylim=ylim, n_levels=256)
    
    # get the color code for red, green and blue
#    red = [graphics.JET256['red'][q] for q in qy]
#    green = [graphics.JET256['green'][q] for q in qy]
#    blue = [graphics.JET256['blue'][q] for q in qy]

#    red = [graphics.COOL256['red'][q] for q in qy]
#    green = [graphics.COOL256['green'][q] for q in qy]
#    blue = [graphics.COOL256['blue'][q] for q in qy]
    

    # rscript
    rscript=''
    if vertical:
        # draw the outer frame
        if x[0]>fromto[0]:
            rscript+='polygon(x=c(%s, %s, %s, %s),y=c(%s, %s, %s, %s))\n' %(str(bot), str(top), str(top), str(bot), str(fromto[0]), str(fromto[0]), str(x[0]), str(x[0]))
#        med=sum(ylim)/2.0
        
        # in the R script, y represents x in our data set because we are drawing a vertical bar
        #rscript += 'y <- c(%s)\n' %str(list(x))[1:-1]
        if x[-1]<fromto[1]:
            rscript+='y<-c(%s, %s)\n' %(str(list(x))[1:-1], str(fromto[1]))
        else:
#            rscript+='y<-c(%s, %s)\n' %(str(list(x))[1:-1], str(x[-1]))
            rscript += 'y<-c(%s)\n' %str(list(x))[1:-1]
        
        
#        rscript += 'cols <- rgb(c(%s), c(%s), c(%s))\n' %(str(red)[1:-1], str(green)[1:-1], str(blue)[1:-1])
        rscript += 'vals <- c(%s)\n' %str(list(y))[1:-1]
        rscript += 'vals[vals > %s] <- %s\n' %(str(ylim[1]), str(ylim[1]))
        rscript += 'vals[vals < %s] <- %s\n' %(str(ylim[0]), str(ylim[0]))
        rscript += 'vals <- round((length(cmap)-1) * (vals - %s)/(%s - %s)) + 1\n' %(str(ylim[0]), str(ylim[1]), str(ylim[0]))
        rscript += 'cols <- cmap[vals]\n'
        rscript += 'for (i in 1:length(cols)) {\n'
        rscript += '\tpolygon(x=c(%s, %s, %s, %s), y=c(y[i], y[i], y[i+1], y[i+1]), col=cols[i], border=cols[i])\n' %(str(bot), str(top), str(top), str(bot))
        rscript += '}\n'
        
#        rscript+='x<-c(%s)\n' %(str(list(y))[1:-1])           
#        rscript+='for (i in 1:length(x)) {\n'
#        rscript+='\tif (x[i]<=%f) {\n' %med
#        rscript+='\t\tcol<-rgb(0, log10((%f-max(x[i],%f))/(%f-%f)*9+1), 0)\n' %(med, ylim[0],ylim[1], ylim[0])
#        rscript+='\t} else {\n'
#        rscript+='\t\tcol<- rgb(log10((min(x[i],%f)-%f)/(%f-%f)*9+1), 0, 0)\n' %(ylim[1],med, ylim[1], ylim[0])
#        rscript+='\t}\n'
#        rscript+='\tpolygon(x=c(%d,%d,%d,%d),y=c(y[i],y[i],y[i+1],y[i+1]),col=col,border=col)\n' %(bot,top,top,bot)
#        rscript+='}\n'
    else:
        if x[0]>fromto[0]:
            rscript+='polygon(x=c(%s, %s, %s, %s), y=c(%s, %s, %s, %s))\n' %(str(fromto[0]), str(x[0]), str(x[0]), str(fromto[0]), str(bot), str(bot), str(top), str(top))
    
#        med=sum(ylim)/2.0
    
        if x[-1]<fromto[1]:
            rscript+='x<-c(%s, %s)\n' %(str(list(x))[1:-1], str(fromto[1]))
        else:
#            rscript+='x<-c(%s, %s)\n' %(str(list(x))[1:-1], str(x[-1]))
            rscript += 'x<-c(%s)\n' %str(list(x))[1:-1]
        
        # quantize the values and get the color code
        #rscript += 'cols <- rgb(c(%s), c(%s), c(%s))\n' %(str(red)[1:-1], str(green)[1:-1], str(blue)[1:-1])
        rscript += 'vals <- c(%s)\n' %str(list(y))[1:-1]
        rscript += 'vals[vals > %s] <- %s\n' %(str(ylim[1]), str(ylim[1]))
        rscript += 'vals[vals < %s] <- %s\n' %(str(ylim[0]), str(ylim[0]))
        rscript += 'vals <- round((length(cmap)-1) * (vals - %s)/(%s - %s)) + 1\n' %(str(ylim[0]), str(ylim[1]), str(ylim[0]))
        rscript += 'cols <- cmap[vals]\n'
        rscript += 'for (i in 1:length(cols)) {\n'
        rscript += '\tpolygon(x=c(x[i], x[i+1],x [i+1], x[i]),y=c(%s, %s, %s, %s), col=cols[i], border=cols[i])\n' %(str(bot), str(bot), str(top), str(top))
        rscript += '}\n'
        
#        rscript+='y<-c(%s)\n' %str(list(y))[1:-1]
#        rscript+='for (i in 1:length(y)) {\n'
#        rscript+='\tif (y[i]<=%f) {\n' %med
#        rscript+='\t\tcol<- rgb(0, log10((%f-max(y[i],%f))/(%f-%f)*9+1), 0)\n' %(med, ylim[0],ylim[1], ylim[0])
#        rscript+='\t} else {\n'
#        rscript+='\t\tcol<- rgb(log10((min(y[i],%f)-%f)/(%f-%f)*9+1), 0, 0)\n' %(ylim[1],med, ylim[1], ylim[0])
#        rscript+='\t}\n'
#        rscript+='\tpolygon(x=c(x[i],x[i+1],x[i+1],x[i]),y=c(%d,%d,%d,%d),col=col,border=col)\n' %(bot,bot,top,top)
#        rscript+='}\n'

    return rscript


def heatmap_rectangles(xstart, xend, y, xlim, ylim, bot=0, top=1, xaxt = "n", cmap="cmap"):
    """Draw a bar of heatmap.
    
    Parameters;
    1. xstart: the start positions of rectangles
    2. xend: the end positions of rectangles
    3. y: the y values (scores) of rectangles
    3. xlim: x limits
    4. ylim: y limits
    5. bot: bottom of the bar
    6. top: top of the bar
    7. cmap: the colormap variable (vector) name in the R script. WARNING: to use heatmap_bar, a colormap must be given prior to this function and correctly referred to. 
    """
 
    # error checking: check whether xstart, xend and y are iterable types
    rscript = ''
    type_xstart = type(xstart)
    type_xend = type(xend)
    type_y = type(y)
    if type_xstart != list and type_xstart != tuple and type_xstart != array:
        return rscript
    if type_xend != list and type_xend != tuple and type_xend != array:
        return rscript
    if type_y != list and type_y != tuple and type_y != array:
        return rscript
    
    rscript += plot(xlim, [bot,top], xlim=xlim, ylim=[bot,top], tp="n", frame=False, yaxt="n", xaxt=xaxt, xlab="", ylab="", main="")
    # copy the rectangle coordinates and their values. The values are saturated wrt ylim
    if type_xstart == array:
        rscript += 'start <- c(%s)\n' %str(list(xstart))[1:-1]
    else:
        rscript += 'start <- c(%s)\n' %str(xstart)[1:-1]
    if type_xend == array:  
        rscript += 'end <- c(%s)\n' %str(list(xend))[1:-1]
    else:
        rscript += 'end <- c(%s)\n' %str(xend)[1:-1]
    if type_y == array:
        rscript += 'vals <- c(%s)\n' %str(list(y))[1:-1]
    else:
        rscript += 'vals <- c(%s)\n' %str(y)[1:-1]
    rscript += 'vals[vals > %s] <- %s\n' %(str(ylim[1]), str(ylim[1]))
    rscript += 'vals[vals < %s] <- %s\n' %(str(ylim[0]), str(ylim[0]))
    rscript += 'vals <- round((length(cmap)-1) * (vals - %s)/(%s - %s)) + 1\n' %(str(ylim[0]), str(ylim[1]), str(ylim[0]))
    rscript += 'cols <- %s[vals]\n' %cmap
    rscript += 'for (i in 1:length(cols)) {\n'
    rscript += '\tpolygon(x=c(start[i], end[i], end[i], start[i]), y=c(%d, %d, %d, %d), col=cols[i], border=cols[i])\n' %(bot, bot, top, top)
    rscript += '}\n'

    return rscript


def rectangles_with_heights(xstart, xend, y, ylim, bot, top, xaxt = "n", col = ["red"]):
    """Draw a bar of heatmap.
    
    Parameters;
    1. xstart: the start positions of rectangles
    2. xend: the end positions of rectangles
    3. y: the y values (scores) of rectangles
    4. ylim: y limits
    5. bot: bottom of the bar
    6. top: top of the bar
    7. xaxt: "n" = not drawing the x axis; other values = drawing
    8. col: the color of rectangles 
    """
    
    
    # error checking: check whether xstart, xend and y are iterable types
    rscript = ''
    type_xstart = type(xstart)
    type_xend = type(xend)
    type_y = type(y)
    if type_xstart != list and type_xstart != tuple and type_xstart != array:
        return rscript
    if type_xend != list and type_xend != tuple and type_xend != array:
        return rscript
    if type_y != list and type_y != tuple and type_y != array:
        return rscript
    
    # to distinguish if a vaiable name is given for 'col'
    if type(col) == tuple or type(col) == list:
        colstr = 'c(' + ','.join(map(lambda s: '"%s"' %str(s), col)) + ')'
    elif type(col) == str:
        colstr = col
        
    # med point
    med = 1.0*(bot + top)/2
    hdelta = 1.0*(top-bot)/2
    # copy the rectangle coordinates and their values. The values are saturated wrt ylim
    if type_xstart == array:
        rscript += 'start <- c(%s)\n' %str(list(xstart))[1:-1]
    else:
        rscript += 'start <- c(%s)\n' %str(xstart)[1:-1]
    if type_xend == array:  
        rscript += 'end <- c(%s)\n' %str(list(xend))[1:-1]
    else:
        rscript += 'end <- c(%s)\n' %str(xend)[1:-1]
    if type_y == array:
        rscript += 'vals <- c(%s)\n' %str(list(y))[1:-1]
    else:
        rscript += 'vals <- c(%s)\n' %str(y)[1:-1]
    
    # saturate the values
    rscript += 'vals[vals > %s] <- %s\n' %(str(ylim[1]), str(ylim[1]))
    rscript += 'vals[vals < %s] <- %s\n' %(str(ylim[0]), str(ylim[0]))
    
    # normalize the values
    if ylim[0] < 0: 
        rscript += 'heights <- %f * ((vals - %s)/(%s - %s) -0.5) + %s\n' %(hdelta/0.5, str(ylim[0]), str(ylim[1]), str(ylim[0]), str(med))
        rscript += 'for (i in 1:length(heights)) {\n'
        rscript += '\tpolygon(x=c(start[i], end[i], end[i], start[i]), y=c(%s, %s, heights[i], heights[i]), col=%s, border=%s)\n' %(str(med), str(med), colstr, colstr)
    else:
        rscript += 'heights <- %f * ((vals - %s)/(%s - %s)) + %s\n' %(2*hdelta, str(ylim[0]), str(ylim[1]), str(ylim[0]), str(bot))
        rscript += 'for (i in 1:length(heights)) {\n'
        rscript += '\tpolygon(x=c(start[i], end[i], end[i], start[i]), y=c(%s, %s, heights[i], heights[i]), col=%s, border=%s)\n' %(str(bot), str(bot), colstr, colstr)
    rscript += '}\n'
    
    return rscript


def rectangles_with_heights_and_colors(xstart, xend, y, ylim, bot, top, xaxt = "n", cmap='cmap'):
    """Draw bars with heights and colors defined by cmap.
    
    Parameters;
    1. xstart: the start positions of rectangles
    2. xend: the end positions of rectangles
    3. y: the y values (scores) of rectangles
    4. ylim: y limits
    5. bot: bottom of the bar
    6. top: top of the bar
    7. cmap: the colormap variable (vector) name in the R script. WARNING: to use heatmap_bar, a colormap must be given prior to this function and correctly referred to. 
    """
 
    # error checking: check whether xstart, xend and y are iterable types
    rscript = ''
    type_xstart = type(xstart)
    type_xend = type(xend)
    type_y = type(y)
    if type_xstart != list and type_xstart != tuple and type_xstart != array:
        return rscript
    if type_xend != list and type_xend != tuple and type_xend != array:
        return rscript
    if type_y != list and type_y != tuple and type_y != array:
        return rscript
    
    # med point
    med = 1.0*(bot + top)/2
    hdelta = 1.0*(top-bot)/2
    # copy the rectangle coordinates and their values. The values are saturated wrt ylim
    if type_xstart == array:
        rscript += 'start <- c(%s)\n' %str(list(xstart))[1:-1]
    else:
        rscript += 'start <- c(%s)\n' %str(xstart)[1:-1]
    if type_xend == array:  
        rscript += 'end <- c(%s)\n' %str(list(xend))[1:-1]
    else:
        rscript += 'end <- c(%s)\n' %str(xend)[1:-1]
    if type_y == array:
        rscript += 'vals <- c(%s)\n' %str(list(y))[1:-1]
    else:
        rscript += 'vals <- c(%s)\n' %str(y)[1:-1]
    
    # saturate the values
    rscript += 'vals[vals > %s] <- %s\n' %(str(ylim[1]), str(ylim[1]))
    rscript += 'vals[vals < %s] <- %s\n' %(str(ylim[0]), str(ylim[0]))
    rscript += 'quantized.vals <- seq(%s, %s, by=(%s-%s)/(length(%s)-1))\n' %(str(ylim[0]), str(ylim[1]), str(ylim[1]), str(ylim[0]), cmap)
    
    # normalize the values
   
    rscript += 'ixs <- round((length(%s)-1) * (vals - %s)/(%s - %s)) + 1\n' %(cmap, str(ylim[0]), str(ylim[1]), str(ylim[0]))
    rscript += 'for (i in 1:length(vals)) {\n'
    if ylim[0] < 0:
        rscript += '\tjs <- (length(quantized.vals)/2):ixs[i]\n'
        rscript += '\tprev <- %f * ((quantized.vals[js[1]] - %s)/(%s - %s) -0.5) + %s\n' %(hdelta/0.5, str(ylim[0]), str(ylim[1]), str(ylim[0]), str(med))
        rscript += '\tfor(j in js) {\n'
        rscript += '\t\theight <- %f * ((quantized.vals[j] - %s)/(%s - %s) -0.5) + %s\n' %(hdelta/0.5, str(ylim[0]), str(ylim[1]), str(ylim[0]), str(med))
        rscript += '\t\tpolygon(x=c(start[i], end[i], end[i], start[i]), y=c(prev, prev, height, height), col=%s[j], border=%s[j])\n' %(cmap, cmap)
        rscript += '\t\tprev <- height\n'
    else:
        rscript += '\tjs <- 1:ixs[i]\n'
        rscript += '\tprev <- %f * ((quantized.vals[js[1]] - %s)/(%s - %s)) + %s\n' %(2*hdelta, str(ylim[0]), str(ylim[1]), str(ylim[0]), str(bot))
        rscript += '\tfor (j in js) {\n'
        rscript += '\t\theight <- %f * ((quantized.vals[j] - %s)/(%s - %s)) + %s\n' %(2*hdelta, str(ylim[0]), str(ylim[1]), str(ylim[0]), str(bot))
        rscript += '\t\tpolygon(x=c(start[i], end[i], end[i], start[i]), y=c(prev, prev, height, height), col=%s[j], border=%s[j])\n' %(cmap, cmap)
        rscript += '\t\tprev <- height\n'
    rscript += '\t}\n'
    rscript += '}\n'
    
    return rscript


def write_func_quantize():
    """Write a quantize function as a R script.
    
    To use quantize function in R, use use_func_quantize.
    """
    
    rscript = 'quantize <- function(x, lim, n.levels)\n'
    rscript += '\t{\n'
    rscript += '\t\t' + 'x[x < lim[0]] <- lim[0]' + '\n'
    rscript += '\t\t' + 'x[x > lim[1]] <- lim[1]' + '\n'
    rscript += '\t\t' + 'q <- round(n.levels * (x - lim[0])/(lim[1]-lim[0]))' + '\n'
    rscript += '\t\t' + 'return(q)' + '\n'
    rscript += '\t\}\n'
    
    return rscript

def use_func_quantize(x, lim, n_levels='128', return_varname = 'qx'):
    """ Quantize x with the given number of levels
    
    parameters:
    1. x: the name of vector or string of R vector
    2. lim: the name of a vector of lower and upper limits or R vector
    3. n_levels: the name of variable of the number of levels or a string number
    4. return_varname: the return argument name in the R script
    """
    
    rscript = '%s <- quantize(%s, lim=%s, n.levels=%s)\n' %(return_varname, x, lim, n_levels)
    
    return rscript
