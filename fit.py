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
from datetime import date
import defaults
import lsqfit, copy, math
import yaml
from myclasses import Ensemble, Fit2
d = date.today()

pfile = open(sys.argv[1],'r')
rfile = open(sys.argv[2],'r')

def main(pfile,rfile):
    # Get ensemble parameters and run parameters from paramfiles
    par = Ensemble(pfile)
    run = yaml.load(rfile)
    no_fits = run['fits']

    #print(run)
    for ifit in range(no_fits):
        fit = run[str(ifit)]
        dofit(par,run,fit)
###

def dofit(params,run,thisfit):

    if thisfit['type'] == '2pt':
        fit = Fit2(thisfit,params)
    fit.results = fit.fitter.lsqfit(data=params.data,prior=fit.prior,p0=None,print_fit=False)
    fit.print_fit(params,'results/'+fit.name+'.'+d.isoformat()+'.fit')
    if fit.nexp == 4:
        fit.print_prior(params,'results/'+fit.name+'.'+d.isoformat()+'.prior')
        fit.print_model(params,'results/'+fit.name+'.'+d.isoformat()+'.model')
        chi2 = fit.chi2_cal()
        chi2_aug = fit.chi2_aug_part()
###

if __name__ == '__main__':
    if True:
        main(pfile,rfile)
