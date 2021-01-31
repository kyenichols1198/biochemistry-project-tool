# -----------------------------------------------------------
# demonstrates how to model the original matlab code
# email kdnichols@gmail.com
# This project was created for a class at UW-Madison
# -----------------------------------------------------------

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import brentq


YSCALE=.0028
def file_read(fname):
        content_array = []
        with open(fname) as f:
                #Content_list is the list that contains the read lines.
                for line in f:
                        content_array.append(YSCALE*float(line))
                return content_array
def convertAraToY (induction):
    x=induction
    # c=333.8
    c=0
    k=505.7
    a=4286
    n=1.775
    return a*(x**n)/(k+(x**n))+c



def hysteresis_test():
    paramBeta = 1773.574498
    paramK = 2952.672305
    araConv=0.036299

    SHIFT=1
    range_dom=[]
    range_exp=[-2,2.4]
    step=.001
    i = range_exp[0]
    while i < range_exp[1]:
        range_dom.append(10**i)
        i+=step
        
    lowest=0
    highest=len(range_dom)-1
    topBranch=[]
    topBranch_dom=[]
    alpha = 10**-5
    def myfun(x,y):
        return alpha+paramBeta*(x-1-y+np.sqrt((x+y+1)**2-4*x*y))/(paramK+x-1-y+np.sqrt((x+y+1)**2-4*x*y))-x
    oldx=fsolve(myfun,1000, args=convertAraToY(lowest)*araConv)[0]

    sn1=len(range_dom)-1
    for idx in range(0,len(range_dom)-1):
        y=convertAraToY(range_dom[idx])*araConv
        x=fsolve(myfun,oldx,args=y)[0]
        #if (x-oldx)**2 < 5*x:
        if idx != 2953:
            topBranch.append(x)
            topBranch_dom.append(range_dom[idx])
            oldx=x
        else:
            topBranch.append(x)
            topBranch_dom.append(range_dom[idx])
            sn1=idx
            break

    sn2=0
    y = convertAraToY(highest)*araConv
    oldx=fsolve(myfun,0,args=y)[0]

    if np.isfinite(oldx):
        bottomBranch=[]
        bottomBranch_dom=[]
        for idx in reversed(range(0,len(range_dom)-1)):
            y=convertAraToY(range_dom[idx])*araConv
            x=fsolve(myfun,oldx,args=y)[0]
            #if (x-oldx)**2 < x/10:
            if idx != 1903:
                bottomBranch.append(x)
                bottomBranch_dom.append(range_dom[idx])
                
                oldx=x
            else:
                sn2=idx
                break

    unstableUpperBound=topBranch[sn1]
    unstableLowerBound=bottomBranch[sn2]
    unstableBranch=[]
    unstableBranch_dom=[]

    y = convertAraToY(range_dom[1903])*araConv
    oldx=fsolve(myfun,bottomBranch[len(bottomBranch)-1],args=y)[0]

    for idx in range(1903,2950):
        y=convertAraToY(range_dom[idx])*araConv
        x= brentq(myfun,oldx,unstableUpperBound-.001,args=y)
        unstableBranch.append(x)
        unstableBranch_dom.append(range_dom[idx])
        oldx = x

    plt.plot(np.array(topBranch_dom), np.array(topBranch),label="upper stable")
    plt.plot(np.array(bottomBranch_dom), np.array(bottomBranch), label="unstable")
    plt.plot(np.array(unstableBranch_dom), np.array(unstableBranch), label="lower stable")

    plt.axvline(x=10**-1, label="monostable 1")
    plt.axvline(x=10**0.5, label="bistable")
    plt.axvline(x=10**1.5, label="monostable 2")

    plt.title('hysteresis diagram')
    plt.xlabel('Total anti-sigma factor')
    plt.xscale('log')
    plt.ylabel('total sigma factor')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.show()



