#!/user/bin/env python
# encoding: utf-8
"""
myexample.py  --- my first attempt to use this package to do a simple 2pt corr fn fit

Written by Thomas Primer 2013-9-11 mostly copied from corr2corr3.py by Peter Lepage
"""

from __future__ import print_function # makes it work for python2 and 3 apparently

import os, sys
import json
import collections
from fitutils import my_3pt_print_results, my_build_3pt_prior, save_results_3pt
from corrfitter import Corr2, CorrFitter, fastfit, Corr3
from gvar import gvar, log, exp, BufferDict, fmt_errorbudget
from gvar.dataset import Dataset, avg_data, bin_data
from numpy import array, arange
import lsqfit
import copy
import math

svdcut = (1e-4,1e-4)

DISPLAYPLOTS = False
try:
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False # in case it doesn't work

stub = sys.argv[1]
2pt1 = 
global tmin
global tmax

def main():

    global tmin
    global tmax

    ifile = open("inputs/"+stub+".input",'r')
    dfile = ifile.readline()
    nt = int(ifile.readline())
    tmin = int(ifile.readline())
    tmax = int(ifile.readline())
    PIMAX = int(ifile.readline())
    DMAX = int(ifile.readline())
    bin = int(ifile.readline())
    exp = int(ifile.readline())
    oexp = int(ifile.readline())
    ncon = int(ifile.readline())
    ntext = int(ifile.readline())
    jack = int(ifile.readline())
    T1 = int(ifile.readline())
    T2 = int(ifile.readline())
    T3 = int(ifile.readline())
    T4 = int(ifile.readline())
    T5 = int(ifile.readline())
    ml = float(ifile.readline())
    ms = float(ifile.readline())
    mc = float(ifile.readline())    
    pi0 = float(ifile.readline())
    pi0_err = float(ifile.readline())
    mpi = gvar(pi0,pi0_err)

    if jack >0:
        print("Doing jackknife\n")
        data = avg_data( bin_data(Dataset(dfile.rstrip()),bin) ,jack=jack-1)
    else:
        data_full = Dataset(dfile.rstrip())
        data = avg_data(bin_data(data_full,bin))
    p0pi=None
    p0D=None
    #p0file="p0.DtoPI"
    savep0=None
    p0file=None
    for nexp in range(1,exp+1):
        if nexp == 1:
            ostart = 1
        else:
            ostart = nexp - 1
        for noexp in range(ostart,nexp+1):
            # Fit first particle
            pifitter = CorrFitter(models=build_models_2pt(2pt1,tmin,tmax,nt,nexp,noexp))
            prior = my_build_prior(stub, nexp, noexp, 'pi')
            pifit = pifitter.lsqfit(data=data, prior=prior, p0=p0pi, savep0file=savep0, ncon=ncon)
            my_print_results(pifit, prior, data, exp, oexp, nexp, noexp, tmin, tmax, bin, 2pt1 )
            save_priors(fit,nexp,noexp,2pt1+"this.prior")
            p0pi = fit.pmean
            # Fit second particle
            prior = my_build_prior(stub, nexp, noexp, 'D')
            fit = fitter.lsqfit(data=data, prior=prior, p0=p0D, savep0file=savep0, ncon=ncon)
            save_priors(fit,nexp,noexp,2pt1+"this.prior"
            p0D = fit.pmean
    # Fit 3pt
    fitter = CorrFitter(models=build_models(tmin,tmax,PIMAX,DMAX,nt,ntext,T1,T2,T3,T4,T5))
    my_3pt_print_results(fit, prior, data, exp, exp, nexp, noexp, tmin, tmax, bin, ml, mc, pi0, pi0_err, 'Pion2', 'Dmeson')
    print('Saving results 3pt')
            save_results_3pt("DtoPI",fit,jack,nexp,noexp,bin,ml,mc,pi0,pi0_err)
            if DISPLAYPLOTS:
                fitter.display_plots()

    fitter.print_fitdata(tmin,tmax)
    keys = fitter.keys
    for k in keys:
        for tt in range(0,len(data[k])):
                print('||  {:10s} | {:3d} {:14.7e} {:14.7e}'.format(k,tt,data[k][tt].mean,data[k][tt].sdev))

##

def build_models_2pt(2pt,tmin,tmax,nt,nexp,noexp):
    """ build models """
    tdata = range(nt)
    tp=nt

    if noexp > 0:
        pia=('pi:a','pi:ao')
        pidE=('pi:dE','pi:dEo')
        Von='Von'
        Voo='Voo'
    else:
        pia=('pi:a')
        pidE=('pi:dE')

    2pt = Corr2( datatag='2pt'+2pt, tp=tp, tdata=tdata, tfit=range(tmin,PIMAX),
                  a=pia, b=pia,dE=pidE, s=(1.,-1.) )

    models = [ 2pt ]
        
#    return [models[:2]] + models[2:]
    return models
##
def build_models_3pt(2pt1,2pt2,3pt,tmin,tmax,PIMAX,DMAX,nt,nexp,noexp,ntext,T1,T2,T3,T4,T5):
    """ build models """
    tdata = range(nt)
    tp=nt

    if noexp > 0:
        pia=('pi:a','pi:ao')
        pidE=('pi:dE','pi:dEo')
        Von='Von'
        Voo='Voo'
    else:
        pia=('pi:a')
        pidE=('pi:dE')
        Von=None
        Voo=None

    Da=('D:a','D:ao')
    DdE=('D:dE','D:dEo')
    Vnn='Vnn'
    Vno='Vno'

    2pt1 = Corr2( datatag='2pt'+2pt1, tp=tp, tdata=tdata, tfit=range(tmin,PIMAX),
                  a=pia, b=pia,dE=pidE, s=(1.,-1.) )

    2pt2 = Corr2( datatag='2pt'+2pt2, tp=tp, tdata=tdata, tfit=range(tmin,DMAX),
                  a=Da, b=Da,dE=DdE, s=(1.,-1.) )

    3ptT1 = Corr3( datatag='3pt'+3pt+str(T1),tdata=range(T1+1),T=T1,tfit=range(tmin,T1-tmax),
                  a=pia, dEa=pidE, sa=(1.,-1), b=Da, dEb=DdE, sb=(1.,-1),
                  Vnn=Vnn,Vno=Vno,Von=Von,Voo=Voo )

    3ptT2 = Corr3( datatag='3pt'+3pt+str(T2),tdata=range(T2+1),T=T2,tfit=range(tmin,T2-tmax),
                  a=pia, dEa=pidE, sa=(1.,-1), b=Da, dEb=DdE, sb=(1.,-1),
                  Vnn=Vnn,Vno=Vno,Von=Von,Voo=Voo )

    3ptT3 = Corr3( datatag='3pt'+3pt+str(T3),tdata=range(T3+1),T=T3,tfit=range(tmin,T3-tmax),
                  a=pia, dEa=pidE, sa=(1.,-1), b=Da, dEb=DdE, sb=(1.,-1),
                  Vnn=Vnn,Vno=Vno,Von=Von,Voo=Voo )

    3ptT4 = Corr3( datatag='3pt'+3pt+str(T4),tdata=range(T4+1),T=T4,tfit=range(tmin,T4-tmax),
                  a=pia, dEa=pidE, sa=(1.,-1), b=Da, dEb=DdE, sb=(1.,-1),
                  Vnn=Vnn,Vno=Vno,Von=Von,Voo=Voo )

    3ptT5 = Corr3( datatag='3pt'+3pt+str(T5),tdata=range(T5+1),T=T5,tfit=range(tmin,T5-tmax),
                  a=pia, dEa=pidE, sa=(1.,-1), b=Da, dEb=DdE, sb=(1.,-1),
                  Vnn=Vnn,Vno=Vno,Von=Von,Voo=Voo )

    if(ntext == 5):
        models = [ 2pt1, 2pt2, 3ptT1, 3ptT2, 3ptT3, 3ptT4, 3ptT5 ]
    if(ntext == 4):
        models = [ 2pt1, 2pt2, 3ptT1, 3ptT2, 3ptT3, 3ptT4 ]
    if(ntext == 3):
        models = [ 2pt1, 2pt2, 3ptT1, 3ptT2, 3ptT3 ]
    if(ntext == 2):
        models = [ 2pt1, 2pt2, 3ptT1, 3ptT2 ]
    if(ntext == 1):
        models = [ 2pt1, 2pt2, 3ptT1 ]
        
#    return [models[:2]] + models[2:]
    return models
##

def fmtlist(x):
    """ Make neat lists for printing. """
    return '  '.join([xi.fmt() for xi in x])
##


def print_ffit(meson,  ffit,  fit):
    osc = ffit.osc
    elabel = "Eo" if osc else "E"
    alabel = "ao" if osc else "a"
    E = ffit.E 
    if osc:
        E_fit = exp(fit.p["log(%s:dEo)" % meson][0])
        a_fit = exp(fit.p["log(%s:ao)" % meson][0])
    else:
        E_fit = exp(fit.p["log(%s:dE)" % meson][0])
        a_fit = exp(fit.p["log(%s:a)" % meson][0])
        a = (ffit.ampl)**0.5
        chi2_dof = ffit.chi2/ffit.dof
        Q = ffit.Q
        print("%8s: %2s = %s   %2s_fit = %s   chi2/dof = %.2f   Q = %.1f"
              % (meson,  elabel,  E.fmt(),  elabel,  E_fit.fmt(),  chi2_dof[0],  Q[0]))
        if osc:
            print("")
        else:
            print("%8s  %2s = %s   %2s_fit = %s   chi2/dof = %.2f   Q = %.1f\n"
                  % ("",  alabel,  a.fmt(),  alabel,  a_fit.fmt(),  chi2_dof[1],  Q[1]))
##

if __name__ == '__main__':
    if True:
        main()
    else:
        import hotshot, hotshot.stats
        prof = hotshot.Profile("lsqfit-test1.prof")
        prof.runcall(main)
        prof.close()
        stats = hotshot.stats.load("lsqfit-test1.prof")
        stats.strip_dirs()
        stats.sort_stats('time', 'calls')
        stats.print_stats(140)
