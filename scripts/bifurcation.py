# -----------------------------------------------------------
# demonstrates how to model the original matlab code
# email kdnichols@gmail.com
# This project was created for a class at UW-Madison
# -----------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def find_saddle(f,b,k,alpha):
    x,y = f
    free=(x-1-y)+np.sqrt((x+y+1)**2-4*x*y)
    der=1+(-4*y+2*(1+x+y))/(2*np.sqrt(-4*x*y+(1+x+y)**2))
    deriv=b*der/(k+free)-b*der*free/((k+free)**2)
    dx=alpha+b*free/(k+free)-x
    return [dx, deriv-1]
    
def Y_to_Ara(y):
    c=0
    k=505.7;
    a=4286;
    n=1.775;
    h=n;
    if y>a:
        ara=10**6;
    else:
        ara=(-1*(y-a)/(k*y))**(-1/h);

    return ara

def aTc_to_alpha(atc,scale):
    x=atc
    a=1990
    c=18.02
    c=1.7
    k=.2286
    h=.5588
    
    alpha=c+a*(x**h)/(k+(x**h))
    return scale*alpha

def in_function(test_vector,vals1,vals2):
    yes=False
    if len(vals1) > 0 and len(vals2) > 0:
        for i in range(0, len(vals1)):
            if (test_vector[0]-vals1[i-1])**2 < 1e-6 and (test_vector[1]-vals2[i-1])**2 < 1e-6:
                yes=True
    return yes
    
def find_roots(b,k,alpha):
    xVals=[]
    yVals=[]
    i =0
    STEPS=[1, .5, .25, .1, .01]
    step =0
    while i <= 4:
        step = STEPS[i]
        while i <= 10:
            soln=fsolve(find_saddle, [10**i, 10**i], args=(b,k,alpha))
            
            if soln[0] >= 0 and soln[1] >=0 and soln[1] < 20 and np.around(soln[0], decimals=2) not in np.around(xVals, decimals=2) and np.around(soln[1], decimals=2) not in np.around(yVals, decimals=2):
                xVals.append(soln[0])
                yVals.append(soln[1])
            i= i + step
        i = i +1
    return [xVals, yVals]

def bifurcation_test():
    beta=1773.574498
    paramk=2952.672305
    araConv=0.036299
    atcConv=.194
    step=.01
    atc_exponent=[-6.5,-2.7]
    atcRange=[]
    i = -6.5

    while i < -2.7:
        atcRange.append(10**i)
        i =i+step

    z = aTc_to_alpha(atcRange[0], atcConv)
    xV = find_roots(beta, paramk, z)[0]
    sn = find_roots(beta, paramk, aTc_to_alpha(atcRange[0], atcConv))[1]

    s1=sn[0]
    s2=sn[1]
    prev1=[xV[0],s1]
    prev2=[xV[1],s2]
    line_r =[]
    line_l =[]

    for i in atcRange:
        alpha = aTc_to_alpha(i, atcConv)
        x1=fsolve(find_saddle,prev1, args=(beta,paramk,alpha))
        x2=fsolve(find_saddle,prev2, args=(beta,paramk,alpha))
        if np.real(x1.all()) and np.real(x2.all()):
            y=[x1[1],x2[1]]
            line_l.append(Y_to_Ara(np.array(y).min()/araConv))
            line_r.append(Y_to_Ara(np.array(y).max()/araConv))
            prev1=x1
            prev2=x2

    plt.plot(np.array(line_l), np.array(atcRange))
    plt.plot(np.array(line_r), np.array(atcRange))
    plt.title('Bifurcation Diagram')
    plt.xlim(0,18)
    plt.ylim(0,.0006)
    plt.xlabel('Ara (10^-5)')
    plt.ylabel('aTc (uM)')
    plt.show()




  



