#!/user/bin/env python
# encoding: utf-8
#from fit import main
from datetime import date
import yaml

do2ptfits = False

def main():
    d = date.today()
    fitname = 'DP'
    nexp = 4
    c_name = 'pimom2'
    sets = { "a06m0008": {"c_tmins":(2,3,4,5,6,7,8),"c_tmaxs":(30,),"p_tmaxs":(36,),"Texts":(31,39,40)},
             "a09m0012": {"c_tmins":(2,3,4,5,6,7,8),"c_tmaxs":(21,),"p_tmaxs":(30,),"Texts":(23,27,32)} }

    sets = { "a06m0008": {"c_tmins":(6,),"c_tmaxs":(30,),"p_tmaxs":(36,),"Texts":(31,39,40)},
             "a09m0012": {"c_tmins":(2,3,4,5,6,7,8),"c_tmaxs":(21,),"p_tmaxs":(30,),"Texts":(23,27,32)} }

    for ens in sets:
        i = 0
        runfile = "runs/"+ens+"/"+fitname+"."+d.isoformat()+".run"
        rfile = open(runfile,'w')
        fit = sets[ens]
        fitpar = (c_name,fit['c_tmins'],fit['c_tmaxs'],nexp,nexp,'2pt','priors/'+ens+'/'+c_name+'.prior')
        i,child = dumpfits_2pt(rfile,fitpar,i)
        fitpar = ('Dmom0',fit['c_tmins'],fit['p_tmaxs'],nexp,nexp,'2pt','priors/'+ens+'/Dmom0.prior')
        i,parent = dumpfits_2pt(rfile,fitpar,i)
        fitpar = ('DP',fit['Texts'],child,parent,True,'priors/DP.prior')
        i = dumpfits_3pt(rfile,fitpar,i)
        rfile.write('fits: '+str(i))
        rfile.close()

def dumpfits_2pt(rfile,fitpar,i):
    child = []
    for tmax in fitpar[2]:
        for tmin in fitpar[1]:
            p0file = None
            savep0file = 'p0file.tmp'
            for nexp in range(1,fitpar[3]+1):
                run = {str(i):
                           {'name': fitpar[0],
                            'tmin': tmin,
                            'tmax': tmax,
                            'nexp': nexp,
                            'noxp': nexp,
                            'type': fitpar[5],
                            'savep0': savep0file,
                            'p0file': p0file,
                            'priorfile': fitpar[6], 
                            'dofit': do2ptfits}}
                p0file = 'p0file.tmp'
                rfile.write(yaml.dump(run))
                i = i+1
            child.append(i-1)
    return i,child

def dumpfits_3pt(rfile,fitpar,i):
    for child in fitpar[2]:
        for parent in fitpar[3]:
            run = {str(i):
                       {'name':fitpar[0],
                        'Text':fitpar[1],
                        'child_index':child,
                        'parent_index':parent,
                        'include2ptdata':fitpar[4],
                        'priorfile':fitpar[5],
                        'type':'3pt' }}
            rfile.write(yaml.dump(run))
            i = i+1
    return i


if __name__ == '__main__':
    if True:
        main()
