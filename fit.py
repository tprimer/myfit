#!/user/bin/env python
# encoding: utf-8
"""
fit.py  ---  general correlator fitting code for 2 and 3 pt fits

Written by Thomas Primer 2016-1-19
"""

from __future__ import print_function # makes it work for python2 and 3 apparently
import os, sys, json, collections
#from fitutils import my_3pt_print_results, my_build_3pt_prior, save_results_3pt
from corrfitter import Corr2, CorrFitter, fastfit, Corr3
from gvar import gvar, log, exp, BufferDict, fmt_errorbudget
from gvar.dataset import Dataset, avg_data, bin_data
from numpy import array, arange
from datetime import date, time, datetime
from time import strftime
import defaults
import lsqfit, copy, math
import yaml
from myclasses import Ensemble, Fit2, Fit3
d = date.today()
t = datetime.now().time()


ens = sys.argv[1]
run_name = sys.argv[2]
#pfile = open(sys.argv[1],'r')
#rfile = open(sys.argv[2],'r')
pfile = open("inputs/"+ens+".params",'r')
rfile = open("runs/"+ens+"/"+run_name+"."+d.isoformat()+".run",'r')

def main(pfile,rfile):
    # Get ensemble parameters and run parameters from paramfiles
    par = Ensemble(pfile)
    run = yaml.load(rfile)
    rfile = 'results/'+ens+'/'+run_name+'.'+d.isoformat()
    no_fits = run['fits']
    fits = []
    p0=None

    for ifit in range(no_fits):
        fit = run[str(ifit)]
        fits.append( dofit(par,run,fit,fits,p0) )
        if fits[ifit].dofit:
            p0 = fits[ifit].results.pmean

    print_header(par,rfile+'.fit')
    print_header(par,rfile+'.good')
    for ifit in range(no_fits):
        if fits[ifit].dofit:
            fits[ifit].print_fit(par,rfile+'.fit')
            if fits[ifit].Q > 0.5 and fits[ifit].f_0.sdev < 0.1*fits[ifit].f_0.mean:
                fits[ifit].print_fit(par,rfile+'.good')
            if fits[ifit].nexp == 3 or fits[ifit].nexp ==4:
                fits[ifit].print_prior(par,rfile+'.prior')
                fits[ifit].print_model(par,rfile+'.model')
###

def dofit(params,run,thisfit,fits,p0):
    if thisfit['type'] == '2pt':
        fit = Fit2(thisfit,params)
    if thisfit['type'] == '3pt':
        fit = Fit3(thisfit,params,fits[thisfit['child_index']],fits[thisfit['parent_index']])
    if fit.dofit:
        fit.results = fit.fitter.lsqfit(data=params.data,prior=fit.prior,p0=None,print_fit=False)
    #fit.print_fit(params,'results/'+fit.name+'.'+d.isoformat()+'.fit')
    return fit
###

def print_header(params,ofile):
    outfile = open(ofile,'a')
    outfile.write("   &    Ranges     & exp &  chi & dof &  q   & Child  &"+
                  " Energy      & Amplit      & Parnt & Energy      & Amplit      & $f_0$        & $q^2$        \\\\ % "+
                  strftime("%H:%M:%S")+" \n")
                  

if __name__ == '__main__':
    if True:
        main(pfile,rfile)
