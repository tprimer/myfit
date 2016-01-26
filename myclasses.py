#!/user/bin/env python
# encoding: utf-8

from corrfitter import Corr2, CorrFitter, fastfit, Corr3
from gvar import gvar, log, exp, BufferDict, fmt_errorbudget, fmt
from gvar.dataset import Dataset, avg_data, bin_data
from pvalue import conf_int
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

        if self.chi2/self.dof > 99:
            chi2dof = 99.9
        else:
            chi2dof = self.chi2/self.dof

        form = "{:s} & {:2d} - {:2d} & {:d}+{:d} & {:4.1f} & {:3d} & {:4.2f} & {:<11s} & {:<11s} \n"
        ofile.write( form.format( self.name,self.tmin,self.tmax,self.nexp,self.noxp,chi2dof,self.dof,self.Q,
                          E[0].fmt(ndecimal=4),a[0].fmt(ndecimal=4) ) )

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
        ofile.write( "   t | model                value                | sigma_m sigma_v \n" )
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
        parameters = len(self.results.p)*self.nexp
        self.dof = len(g)-parameters
        self.chi2 = total
        return total,self.dof

    def chi2_aug_part(self):
        """calculate the part of the chi2 coming from the priors"""
        pth = self.results.p
        p = self.prior
        aug = 0
        for key in pth.keys():
            for exp in range(self.nexp):
                aug = aug + ( (p[key][exp].mean-pth[key][exp].mean)/p[key][exp].sdev )**2
        self.aug = aug
        return aug
        
            

class Fit3:
    
    def __init__(self, dic, childfit, parentfit):
        self.name = dic['name']
        self.Texts = dic['Text']
        self.include2ptdata = dic['include2ptdata']
        self.child = childfit
        self.parent = parentfit
        self.tmin = childfit.tmin
        self.tmax = parentfit.tmin
        self.priorfile = dic['priorfile']
        if childfit.name == 'pimom0':
            Von=None
            Voo=None
        else:
            Von='Von'
            Voo='Voo'
        Vnn='Vnn'
        Vno='Vno'
        self.model = self.build_model(params)
        self.fitter = CorrFitter(self.model)
        self.prior = self.child.prior + self.parent.prior + self.build_prior(params)

            

    def build_3pt_models(self,params):
        models_3pt = []
        for T in Texts:
            models_3pt.append( Corr3( datatag='3pt'+self.name+str(T),tdata=range(T+1),
                                      tfit=range(self.tmin,T-self.tmax),
                                      a=self.child.a,dEa=self.child.dE,sa=(1.,-1.),
                                      b=self.parent.a,dEb=self.parent.dE,sb=(1.,-1.),
                                      Vnn=self.Vnn,Vno=self.Vno,Von=self.Von,Voo=self.Voo) )
        if self.include2ptdata:
            models_3pt.append( child.build_model(params) )
            models_3pt.append( parent.build_model(params) ) 

        return models_3pt
        
    def build_prior(self,params):
        """build interaction matrix prior"""
        prior = BufferDict()
        
        nexp = self.nexp
        noxp = self.noxp
        prior_file = open(self.priorfile,'r')
        pp = yaml.load(prior_file) #
        
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
                    