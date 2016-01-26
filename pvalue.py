""" calculating P-values for correlated fits """

# Written by Thomas Primer
# Date: 2/26/14

import numpy
import gvar
import lsqfit
import time
import collections
import math

def conf_int(chisq, N, DOF):

    #print('***In conf_int***')
    if chisq <= 0.0:
        return 1.0
    if N <= DOF or DOF <= 0:
        return 0.0

    #print('Chisq:',chisq," N:",N," DOF:",DOF)

    if chisq < 3.0*DOF+2.0:
#        print(' R: doing it the first way')
        nsteps = int( 100*chisq/(2.0*(2.0*DOF)**(0.5)) )
        if nsteps < 20:
            nsteps = 20
        eps = chisq/float(nsteps)
        isum=0.0
        chisq_t=0.5*eps
        while ( chisq_t < chisq ):
            isum += Fdist_2( chisq_t, N, DOF )
            chisq_t += eps
        if isum*eps < 0.9:
            return 1.0-isum*eps

        #return 1.0-isum*eps
    

#    print(' R: doing it the second way')
    nsteps = 0
    eps = (2.0*math.sqrt(2.0*DOF))/100.0
    #print(eps)
    y = Fdist_2( chisq, N, DOF )
    
    nsteps = 0
    isum = 0
    chisq_t=chisq+0.5*eps
    while True:
        z = Fdist_2(chisq_t, N, DOF)
        if nsteps > 1000 or z < 0.0001*y:
            break
        isum += z
        nsteps += 1
    return isum*eps

def Fdist_2(chisq, N, DOF):
    
    Nd=float(N)
    Dd=float(DOF)

    a = math.exp( gammln((Nd)/2.0) - gammln(Dd/2.0) - gammln((Nd-Dd)/2.0) \
                 + (-Dd/2.0)*math.log(Nd) \
                 + ((Dd-2.0)/2.0)*math.log(chisq) \
                 + (-Nd/2.0)*math.log(1.0+chisq/Nd) )
    return a

def gammln(xx):
    
    cof = [
        76.18009173,-86.50532033,24.01409822,
        -1.231739516,0.120858003e-2,-0.536382e-5
        ]

    x = xx - 1.0
    tmp = x + 5.5
    tmp -= (x+0.5)*math.log(tmp)
    ser = 1.0
    for j in range(0,6):
        x += 1.0
        ser += cof[j]/x
    return -tmp+math.log(2.50662827465*ser)


