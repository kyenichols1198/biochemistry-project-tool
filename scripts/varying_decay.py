# -----------------------------------------------------------
# demonstrates how to model the original matlab code
# email kdnichols@gmail.com
# This project was created for a class at UW-Madison
# -----------------------------------------------------------

import warnings
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import brentq



def hanging_line(point1, point2, range):


    a = (point2[1] - point1[1])/(np.cosh(point2[0]) - np.cosh(point1[0]))
    b = point1[1] - a*np.cosh(point1[0])
    x = range
    y=[]
    for i in x:
        y.append(a*np.cosh(i) + b)
    return y
def convertAraToY (induction):
    x=induction
    # c=333.8
    c=0
    k=505.7
    a=4286
    n=1.775
    return a*(x**n)/(k+(x**n))+c

def myfun(x,y,paramK, paramBeta, alpha):
    return alpha+paramBeta*(x-1-y+np.sqrt((x+y+1)**2-4*x*y))/(paramK+x-1-y+np.sqrt((x+y+1)**2-4*x*y))-x

def hysteresis_line(paramK, paramBeta, alpha, araConv, gamma):

    SHIFT=1
    YSCALE=.0028
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
    unstableBranch=[]
    unstableBranch_dom=[]
    bottomBranch=[]
    bottomBranch_dom=[]
    
    oldx=fsolve(myfun,1000, args=(convertAraToY(lowest)*araConv,paramK, paramBeta, alpha))[0]

    sn1=len(range_dom)-1
    for idx in range(0,len(range_dom)-1):
        y=convertAraToY(range_dom[idx])*araConv
        x=fsolve(myfun,oldx,args=(y,paramK, paramBeta, alpha))[0]
        #if (x-oldx)**2 < .38*x:
        if idx != 3300:
            topBranch.append(x)
            if topBranch[idx-1] < topBranch[idx]:
                topBranch_dom.append(range_dom[idx])
                sn1=idx
                break
            topBranch_dom.append(range_dom[idx])
            oldx=x



    sn2=0
    y = convertAraToY(highest)*araConv
    oldx=fsolve(myfun,0,args=(y,paramK, paramBeta, alpha))[0]

    if np.isfinite(oldx):
        for idx in reversed(range(0,len(range_dom)-1)):
            y=convertAraToY(range_dom[idx])*araConv
            x=fsolve(myfun,oldx,args=(y,paramK, paramBeta, alpha))[0]
            if (x-oldx)**2 < (x/10):
                bottomBranch.append(x)
                bottomBranch_dom.append(range_dom[idx])
                oldx=x
               
            else:
                bottomBranch.append(x)
                bottomBranch_dom.append(range_dom[idx])
                sn2=idx
                break
    unstableUpperBound=topBranch[sn1]
    unstableLowerBound=bottomBranch[sn2]
    if gamma == 1:

        oldx=fsolve(myfun,unstableLowerBound,args=(y,paramK, paramBeta, alpha))[0]

        for idx in range(sn2,sn1):
            warnings.simplefilter("error")
            try:
                y=convertAraToY(range_dom[idx])*araConv
                x= brentq(myfun,oldx,unstableUpperBound,args=(y,paramK, paramBeta, alpha))
                unstableBranch.append(x)
                unstableBranch_dom.append(range_dom[idx])
                oldx = x
            except (ValueError, RuntimeWarning):
                #unstableUpperBound-=1
                unstableLowerBound+=100
                oldx+=.01
    else:
        unstableBranch = [bottomBranch[len(bottomBranch)-1], topBranch[len(topBranch)-1]]
        unstableBranch_dom = [bottomBranch_dom[len(bottomBranch_dom)-1], topBranch_dom[len(topBranch_dom)-1]]

            


    #for idx in range(1903,2950):
    

    
    #return [topBranch_dom, topBranch, bottomBranch_dom, bottomBranch, x_values, y_values]
    return [topBranch_dom, topBranch, bottomBranch_dom, bottomBranch, unstableBranch_dom,unstableBranch]

################################################################
def tester():
    b = 1773.574498
    gamma_range=[.9,1]
    colors=['red','green','orange','blue','grey','black','yellow']
    color_counter=0
    gamma=gamma_range[0]
    while gamma <= gamma_range[1]:
        

        a =10**-5
        k = 2952.672305
        ara_c=0.036299

        lines =hysteresis_line(k, b*(1/gamma), a, ara_c, gamma)
        line_labels=[]
        print("gamma: "+ str(gamma))
        if gamma !=1:
            upper, = plt.plot(np.array(lines[0]), np.array(lines[1]), color=colors[color_counter], label='gamma= '+str(gamma))
            line_labels.append(upper)
            plt.plot(np.array(lines[2]), np.array(lines[3]), color=colors[color_counter])
            plt.axvline(x=lines[4][0], color=colors[color_counter])
            plt.axvline(x=lines[4][len(lines[4])-1], color=colors[color_counter])
        else:
            upper_original, = plt.plot(np.array(lines[0]), np.array(lines[1]), color=colors[color_counter], label='gamma= '+str(gamma))
            line_labels.append(upper_original)
            plt.plot(np.array(lines[2]), np.array(lines[3]), color=colors[color_counter])
            plt.plot(np.array(lines[4]), np.array(lines[5]), color=colors[color_counter])
            
        gamma+=.02
        color_counter+=1

    ################################################################
    plt.legend()
    plt.title('hysteresis diagram')
    plt.xlabel('Total anti-sigma factor')
    plt.xlim(10**-1,10**2)
    plt.xscale('log')
    plt.ylabel('total sigma factor')
    plt.show()



