#!/user/bin/env python
# encoding: utf-8
#from fit import main
from datetime import date
import yaml

def main():
    d = date.today()
    fitname = '2ptDmom0'
    runfile = "runs/"+fitname+"."+d.isoformat()+".run"
    rfile = open(runfile,'w')
    i=0


    #fitpar = ('pimom2',(2,4,6),(15,),4,4,'2pt','priors/pimom2.prior')
    #i = dumpfits(rfile,fitpar,i)
    fitpar = ('Dmom0',(2,4,6),(20,25,30),4,4,'2pt','priors/Dmom0.prior')
    i = dumpfits(rfile,fitpar,i)

    rfile.write('fits: '+str(i))



def dumpfits(rfile,fitpar,i):
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
                            'priorfile': fitpar[6] }}
                p0file = 'p0file.tmp'
                rfile.write(yaml.dump(run))
                i = i+1
    return i


if __name__ == '__main__':
    if True:
        main()
