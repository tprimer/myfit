#!/user/bin/env python
# encoding: utf-8

from corrfitter import Corr2, CorrFitter, fastfit, Corr3
from gvar import gvar, log, exp, sqrt, BufferDict, fmt_errorbudget, fmt
from gvar.dataset import Dataset, avg_data, bin_data
from pvalue import conf_int
from datetime import date, time, datetime
from time import strftime
import numpy as np
import yaml

class Ensemble:

    def __init__(self, pfile):
        params = yaml.load(pfile)
        self.datafile = params['datafile']
        self.NT = params['NT']
        self.NX = params['NX']
        self.m_l = params['m_l']
        self.m_s = params['m_s']
        self.m_c = params['m_c']
        self.Text = params['Text']
        self.ncon = params['ncon']
        self.pi2twist = params['pi2twist']
        self.k2twist = params['k2twist']
        self.pi0 = params['pi0']
        self.pi0_err = params['pi0_err']
        self.bin = 1
        self.data_full = Dataset(self.datafile.rstrip())
        self.data = avg_data(bin_data(self.data_full,self.bin))

class Fit2:

    def __init__(self, dic, params):
        #print('Creating 2pt fit')
        #print(dic)
        self.name = dic['name']
        self.tmin = dic['tmin']
        self.tmax = dic['tmax']
        self.nexp = dic['nexp']
        self.noxp = dic['noxp']
        self.priorfile = dic['priorfile']
        self.p0file = dic['p0file']
        self.savep0 = dic['savep0']
        self.type = dic['type']
        self.dofit = dic['dofit']
        name = self.name
        if name == 'pimom0':
            self.a=(name+':a')
            self.dE=(name+':dE')
        else:
            self.a=(name+':a',name+':ao')
            self.dE=(name+':dE',name+':dEo')
        # build model
        self.model = self.build_model(params)
        self.fitter = CorrFitter(self.model)
        # build prior
        self.prior = self.build_prior(params)

    def build_model(self,params):
        model = Corr2( datatag='2pt'+self.name,tp=params.NT,tdata=range(params.NT),
                       tfit=range(self.tmin,self.tmax+1),
                       a=self.a,b=self.a,dE=self.dE,s=(1.,-1.) )
        return model

    def build_prior(self,params):
        """build prior from file"""
        prior = BufferDict()

        name = self.name
        nexp = self.nexp
        noxp = self.noxp
        prior_file = open(self.priorfile,'r')
        pp = yaml.load(prior_file)

        prior['log('+name+':dE)'] = [ gvar(0,0) for i in range(nexp) ]
        prior['log('+name+':a)'] = [ gvar(0,0) for i in range(nexp) ]

        # non-osc. state priors
        prior['log('+name+':dE)'][0] = log(gvar(pp['e0'][0],pp['e0'][1]))
        prior['log('+name+':a)'][0] = log(gvar(pp['a0'][0],pp['a0'][1]))
        for i in range(1,nexp):
            prior['log('+name+':dE)'][i] = log(gvar(pp['e1'][0],pp['e1'][1]))
            prior['log('+name+':a)'][i] = log(gvar(pp['a1'][0],pp['a1'][1]))
        # osc. state priors
        if self.name != 'pimom0':
            prior['log('+name+':dEo)'] = [ 0 for i in range(noxp) ]
            prior['log('+name+':ao)'] = [ 0 for i in range(noxp) ]
            prior['log('+name+':dEo)'][0] = log(gvar(pp['o0'][0],pp['o0'][1]))
            prior['log('+name+':ao)'][0] = log(gvar(pp['b0'][0],pp['b0'][1]))
            for i in range(1,noxp):
                prior['log('+name+':dEo)'][i] = log(gvar(pp['o1'][0],pp['o1'][1]))
                prior['log('+name+':ao)'][i] = log(gvar(pp['b1'][0],pp['b1'][1]))
        #print(prior)
        return prior

    def print_fit(self,par,outfile):
        """print the energy and amplitude results from a fit"""
        fit = self.results
        name = self.name
        p = fit.p
        self.chi2_cal()
        self.chi2_aug_part()
        self.Q = conf_int(self.chi2/2, par.ncon, self.dof/2)

        dE = exp(p['log('+name+':dE)'])
        E = [sum(dE[:i+1]) for i in range(self.nexp)]
        a = exp(p['log('+name+':a)'])
        if self.name != 'pimom0':
            dEo = exp(p['log('+name+':dEo)'])
            Eo = [sum(dEo[:i+1]) for i in range(self.nexp)]
            ao = exp(p['log('+name+':ao)'])

        ofile = open(outfile,'a')

        if self.dof > 0:
            chi2dof = self.chi2/self.dof
        else:
            chi2dof = 99.9
        if chi2dof > 99:
            chi2dof = 99.9

        form = "{:s} & {:2d} - {:2d} & {:d}+{:d} & {:4.1f} & {:3d} & {:4.2f} "
        ena  = "& {:<11s} & {:<11s} "
        ofile.write( form.format( self.name,self.tmin,self.tmax,self.nexp,self.noxp,chi2dof,self.dof,self.Q ) )
        for i in range(self.nexp):
            ofile.write( ena.format(E[i].fmt(ndecimal=4),a[i].fmt(ndecimal=4),Eo[i].fmt(ndecimal=4),ao[i].fmt(ndecimal=4)) )
        ofile.write(" \\\\ \n")
                         

    def print_prior(self,par,outfile):
        """print the priors and result for comparison"""
        fit = self.results
        p = self.prior
        name = self.name
        f = fit.p
        
        dE = exp(f['log('+name+':dE)'])
        E = [sum(dE[:i+1]) for i in range(self.nexp)]
        a = exp(f['log('+name+':a)'])
        if self.name != 'pimom0':
            dEo = exp(f['log('+name+':dEo)'])
            Eo = [sum(dEo[:i+1]) for i in range(self.nexp)]
            ao = exp(f['log('+name+':ao)'])
        
        pdE = exp(p['log('+name+':dE)'])
        pa = exp(p['log('+name+':a)'])
        if self.name != 'pimom0':
            pdEo = exp(p['log('+name+':dEo)'])
            pao = exp(p['log('+name+':ao)'])

        ofile = open(outfile,'a')
        ofile.write("{:s} & Prior & Value & & Prior & Value \\\\ # {:d}-{:d} {:d}+{:d} \n".format(
                name,self.tmin,self.tmax,self.nexp,self.nexp))
            
        form = "$a_{:d}$  & {:<11s} & {:<11s} & $dE_{:d}$  & {:<11s} & {:<11s} \\\\ \n"
        for state in range(self.nexp):
            ofile.write( form.format( state,pa[state].fmt(),a[state].fmt(),
                                      state,pdE[state].fmt(),dE[state].fmt() ) )
        if self.name != 'pimom0':
            form = "$a'_{:d}$ & {:<11s} & {:<11s} & $dE'_{:d}$ & {:<11s} & {:<11s} \\\\ \n"
            for state in range(self.nexp):
                ofile.write( form.format( state,pao[state].fmt(),ao[state].fmt(),
                                          state,pdEo[state].fmt(),dEo[state].fmt() ) )
        
    def print_model(self,par,outfile):
        """print the fit data and model and the difference"""
        fit = self.results
        t,g,dg,gth,dgth = self.fitter.collect_fitresults()['2pt'+self.name]
        ofile = open(outfile,'a')
        ofile.write( "   t | {:12s}           value                | sigma_m sigma_v \n".format('2pt'+self.name) )
        for it in range(0,self.tmax-self.tmin):
            data = gvar(g[it],dg[it])
            model = gvar(gth[it],dgth[it])
            diff1 = (data.mean-model.mean)/data.sdev
            diff2 = (data.mean-model.mean)/model.sdev
            ofile.write( " {:3d} | {:<20s} {:<20s} | {:+5.3f}  {:+5.3f} \n".format(
                    t[it],data.fmt(),model.fmt(),diff1,diff2) )

    def chi2_cal(self):
        """calculate chi2 of fit"""
        t,g,dg,gth,dgth = self.fitter.collect_fitresults()['2pt'+self.name]
        chi2 = 0.
        for it in range(0,self.tmax-self.tmin):
            chi2 = ( (g[it]-gth[it])/dg )**2
        total = sum(chi2)
        #print(chi2,total)
        parameters = 0
        for key in self.prior.keys():
            parameters = parameters + len(self.results.p[key])
        self.dof = len(g)-parameters
        self.chi2 = total
        return total,self.dof

    def chi2_aug_part(self):
        """calculate the part of the chi2 coming from the priors"""
        pth = self.results.p
        p = self.prior
        aug = 0
        for key in p.keys():
            for exp in range(self.nexp):
                aug = aug + ( (p[key][exp].mean-pth[key][exp].mean)/p[key][exp].sdev )**2
        self.aug = aug
        return aug
        
            

class Fit3:
    
    def __init__(self, dic, params, childfit, parentfit):
        self.name = dic['name']
        self.Texts = dic['Text']
        self.child_index = dic['child_index']
        self.parent_index = dic['parent_index']
        self.include2ptdata = dic['include2ptdata']
        self.priorfile = dic['priorfile']
        self.child = childfit
        self.parent = parentfit
        self.tmin = childfit.tmin
        self.tmax = parentfit.tmin
        self.nexp = self.child.nexp
        self.dofit = True
        if childfit.name == 'pimom0':
            self.Von=None
            self.Voo=None
        else:
            self.Von='Von'
            self.Voo='Voo'
        self.Vnn='Vnn'
        self.Vno='Vno'
        self.model = self.build_3pt_models(params)
        self.fitter = CorrFitter(self.model)
        self.prior =  self.build_prior(params)
        self.prior.update( self.child.prior )
        self.prior.update( self.parent.prior )

            

    def build_3pt_models(self,params):
        models_3pt = []
        for T in self.Texts:
            models_3pt.append( Corr3( datatag="3pt"+self.name+"T"+str(T),T=T,tdata=range(T+1),
                                      tfit=range(self.tmin,T-self.tmax+1),
                                      a=self.child.a,dEa=self.child.dE,sa=(1.,-1.),
                                      b=self.parent.a,dEb=self.parent.dE,sb=(1.,-1.),
                                      Vnn=self.Vnn,Vno=self.Vno,Von=self.Von,Voo=self.Voo) )
        if self.include2ptdata:
            models_3pt.append( self.child.build_model(params) )
            models_3pt.append( self.parent.build_model(params) ) 

        return models_3pt
        
    def build_prior(self,params):
        """build interaction matrix prior"""
        prior = BufferDict()
        
        nexp = self.child.nexp
        noxp = self.child.noxp
        prior_file = open(self.priorfile,'r')
        pp = yaml.load(prior_file)
        
        prior['Vnn'] = [[ gvar(pp['v22'][0],pp['v22'][1]) for i in range(nexp)] for j in range(nexp)]
        prior['Vno'] = [[ gvar(pp['v22'][0],pp['v22'][1]) for i in range(nexp)] for j in range(nexp)]
        
        prior['Vnn'][0][0] = gvar(pp['v00'][0],pp['v00'][1])
        prior['Vno'][0][0] = gvar(pp['v11'][0],pp['v11'][1])
        for i in range(1,nexp):
            for j in range(1,nexp):
                if i<1 or j<1:
                    prior['Vnn'][i][j] = gvar(pp['v11'][0],pp['v11'][1])
                    prior['Vno'][i][j] = gvar(pp['v11'][0],pp['v11'][1])
                    
        if self.name != "DPm0":
            prior['Voo'] = [[ gvar(pp['v22'][0],pp['v22'][1]) for i in range(nexp)] for j in range(nexp)]
            prior['Von'] = [[ gvar(pp['v22'][0],pp['v22'][1]) for i in range(nexp)] for j in range(nexp)]
        
            prior['Voo'][0][0] = gvar(pp['v11'][0],pp['v11'][1])
            prior['Von'][0][0] = gvar(pp['v11'][0],pp['v11'][1])
            for i in range(1,nexp):
                for j in range(1,nexp):
                    if i<1 or j<1:
                        prior['Voo'][i][j] = gvar(pp['v11'][0],pp['v11'][1])
                        prior['Von'][i][j] = gvar(pp['v11'][0],pp['v11'][1])
                        
        return prior
                    
    def print_fit(self,par,outfile):
        """print the energy and amplitude results from a fit"""
        fit = self.results
        name = self.name
        p = fit.p
        self.chi2_cal()
        #self.chi2_aug_part()
        self.Q = conf_int(self.chi2/2, par.ncon, self.dof/2)
        nexp = self.child.nexp

        name = self.child.name
        CdE = exp(p['log('+name+':dE)'])
        CE = [sum(CdE[:i+1]) for i in range(nexp)]
        a = exp(p['log('+name+':a)'])
        if self.name != 'pimom0':
            CdEo = exp(p['log('+name+':dEo)'])
            CEo = [sum(CdEo[:i+1]) for i in range(nexp)]
            ao = exp(p['log('+name+':ao)'])

        name = self.parent.name
        PdE = exp(p['log('+name+':dE)'])
        PE = [sum(PdE[:i+1]) for i in range(nexp)]
        b = exp(p['log('+name+':a)'])
        if self.name != 'pimom0':
            PdEo = exp(p['log('+name+':dEo)'])
            PEo = [sum(PdEo[:i+1]) for i in range(nexp)]
            bo = exp(p['log('+name+':ao)'])

        # calculating f_0
        # DOESN"T WORK FOR D TO K YET!!!!!!!!!
        if self.child.name == 'pimom0':
            mpi = CE[0]
        else:
            mpi = gvar(par.pi0,par.pi0_err)
        if self.child.name == 'pimom0' or self.child.name == 'pimom2':
            m_q = par.m_l
        else:
            m_q = par.m_s
        Epi = CE[0]
        mD = PE[0]
        v = p['Vnn'][0][0]
        self.f_0 = v*sqrt(Epi*mD)*(par.m_c-m_q)/(mD**2-mpi**2)
        self.qsq = mpi**2+mD**2-2*mD*Epi        

        ofile = open(outfile,'a')

        if self.dof != 0:
            chi2dof = self.chi2/self.dof
        else:
            chi2dof = 99.9
        if chi2dof > 99:
            chi2dof = 99.9

        pars = "{:s} & {:2d}-{:2d} & {:2d}-{:2d} & {:d}+{:d} & {:4.1f} & {:3d} & {:4.2f} & "
        energies = "{:s} & {:<11s} & {:<11s} & "
        form = "{:<12s} & {:<12s} \\\\ \n"
        ofile.write( pars.format( self.name,self.tmin,self.child.tmax,self.tmax,self.parent.tmax,
                                  nexp,nexp,chi2dof,self.dof,self.Q )
                     +energies.format( self.child.name,CE[0].fmt(ndecimal=4),a[0].fmt(ndecimal=4) )
                     +energies.format( self.parent.name,PE[0].fmt(ndecimal=4),b[0].fmt(ndecimal=4) )
                     +form.format( self.f_0.fmt(ndecimal=4), self.qsq.fmt(ndecimal=4) ) )

    def print_model(self,par,outfile):
        """print the fit data and model and the difference"""
        fit = self.results
        names = ["2pt"+self.child.name,"2pt"+self.parent.name]
        [names.append("3pt"+self.name+"T"+str(T)) for T in self.Texts]
        print(names)
        ofile = open(outfile,'a')
        ofile.write("{:5s} models, {:d}+{:d} fit, {:d} to T-{:d} window "+strftime("%H:%M:%S")+"\n".format(self.name,self.nexp,self.nexp,self.child.tmin,self.child.tmax))
        for name in names:
            t,g,dg,gth,dgth = self.fitter.collect_fitresults()[name]
            ofile.write( "   t | {:<11s}          theory               | sigma_data sigma_th \n".format(name) )
            for it in range(0,len(t)):
                data = gvar(g[it],dg[it])
                model = gvar(gth[it],dgth[it])
                diff1 = (data.mean-model.mean)/data.sdev
                diff2 = (data.mean-model.mean)/model.sdev
                ofile.write( " {:3d} | {:<20s} {:<20s} | {:+6.3f}     {:+6.3f} \n".format(
                        t[it],data.fmt(),model.fmt(),diff1,diff2) )

    def print_prior(self,par,outfile):
        """print the priors and result for comparison"""
        fit = self.results
        p = self.prior
        name = self.name
        f = fit.p

        ofile = open(outfile,'a')

        ofile.write("{:5s} priors, {:d}+{:d} fit, {:d} to T-{:d} window "+strftime("%H:%M:%S")+"\n".format(name,self.nexp,self.nexp,self.child.tmin,self.child.tmax))
        self.print_2pt_prior(par,ofile,self.child.name)
        self.print_2pt_prior(par,ofile,self.parent.name)

        ofile.write("V \n")
        form = "{:<8s} "
        self.V_prior_line(ofile,"Vnn prior",form,self.nexp,p['Vnn'])
        self.V_prior_line(ofile,"Vnn fittd",form,self.nexp,f['Vnn'])
        self.V_prior_line(ofile,"Vno prior",form,self.nexp,p['Vno'])
        self.V_prior_line(ofile,"Vno fittd",form,self.nexp,f['Vno'])
        if self.child.name != 'pimom0':
            self.V_prior_line(ofile,"Von prior",form,self.nexp,p['Von'])
            self.V_prior_line(ofile,"Von fittd",form,self.nexp,f['Von'])
            self.V_prior_line(ofile,"Voo prior",form,self.nexp,p['Voo'])
            self.V_prior_line(ofile,"Voo fittd",form,self.nexp,f['Voo'])

    def print_2pt_prior(self,par,ofile,name):
        fit = self.results
        p = self.prior
        f = fit.p
        
        dE = exp(f['log('+name+':dE)'])
        E = [sum(dE[:i+1]) for i in range(self.nexp)]
        a = exp(f['log('+name+':a)'])
        if self.name != 'pimom0':
            dEo = exp(f['log('+name+':dEo)'])
            Eo = [sum(dEo[:i+1]) for i in range(self.nexp)]
            ao = exp(f['log('+name+':ao)'])
        
        pdE = exp(p['log('+name+':dE)'])
        pa = exp(p['log('+name+':a)'])
        if self.name != 'pimom0':
            pdEo = exp(p['log('+name+':dEo)'])
            pao = exp(p['log('+name+':ao)'])

        ofile.write("{:5s} \n".format(name))
        form = "{:<12s} {:<12s} "
        self.prior_line(ofile,"odd prior",form,self.nexp,pa,pdE)
        self.prior_line(ofile,"odd fittd",form,self.nexp,a,dE)
        if self.name != 'pimom0':  
            self.prior_line(ofile,"evn prior",form,self.nexp,pao,pdEo)
            self.prior_line(ofile,"evn fittd",form,self.nexp,ao,dEo)

    def prior_line(self,ofile,name,form,nexp,a,E):
        ofile.write(name+": ")
        ofile.write(form.format(a[0].fmt(),E[0].fmt()))
        for state in range(1,self.nexp):
            ofile.write("| "+form.format(a[state].fmt(),E[state].fmt()))
        ofile.write(" \n")

    def V_prior_line(self,ofile,name,form,nexp,V):
        ofile.write(name+": ")
        ofile.write(form.format(V[0][0].fmt()))
        for i in range(1,self.nexp):
            for j in range(1,self.nexp):
                ofile.write("| "+form.format(V[i][j].fmt()))
        ofile.write(" \n")


    def chi2_cal(self):
        """calculate chi2 of fit"""

        chi2 = 0.
        points = 0
        for fit in (self.child,self.parent):
            t,g,dg,gth,dgth = self.fitter.collect_fitresults()["2pt"+fit.name]
            for it in range(0,fit.tmax-fit.tmin):
                chi2 = ( (g[it]-gth[it])/dg )**2
            total2pt = sum(chi2)
            points = points + len(g)
        for T in self.Texts:
            t,g,dg,gth,dgth = self.fitter.collect_fitresults()["3pt"+self.name+"T"+str(T)]
            for it in range(0,T-fit.tmax-fit.tmin):
                chi2 = ( (g[it]-gth[it])/dg )**2
            total3pt = sum(chi2)
            points = points + len(g)
            
        parameters = 0
        for key in self.prior.keys():
            parameters = parameters + len(self.results.p[key])
        self.dof = points-parameters
        self.chi2 = total2pt+total3pt
        return self.chi2,self.dof

    def chi2_aug_part(self):
        """calculate the part of the chi2 coming from the priors"""
        ####### NOT WORKING ATM #####
        pth = self.results.p
        p = self.prior
        aug = 0
        print(pth.values())
        qth = np.fromiter(pth.values(),np.float)
        q = np.fromiter(p.values(),np.float)
        print(qth)
        print(q)
        for exp in range(self.nexp):
            aug = aug + ( (p[key][exp].mean-pth[key][exp].mean)/p[key][exp].sdev )**2
        self.aug = aug
        return aug
