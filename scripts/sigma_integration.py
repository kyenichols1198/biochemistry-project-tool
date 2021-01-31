# -----------------------------------------------------------
# demonstrates how to model the original matlab code
# email kdnichols@gmail.com
# This project was created for a class at UW-Madison
# -----------------------------------------------------------

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# parameters to switch states
ara = [.1*10e-5,15*10e-5] # metrics are percent concentration
# first position is cell's ON state, and the second OFF state
aTc = .5*120e-2 # uM of transcription factor
# alpha = 0

def convertYtoAra(y):
    c=0  # very "leaky promoter"
    k=505.7
    a=4286
    n=1.775
    h=n
    if y>a:
        ara=10e6
    else:
        ara=(-1*(y-a)/(k*y))**(-1/h)
    return ara

def convertATCtoAlpha(atc,scale):
    
    x=atc
    a=1990
    c=18.02
    k=.2286
    h=.5588
    alpha=c+a*(x**h)/(k+(x**h))
    return scale*alpha
    
def convertAraToY (ind):

    x=ind
    c=0  # very "leaky promoter"
    k=505.7
    a=428
    n=1.775
    y=a*(x**n)/(k+(x**n))+c
    return .0061*y

# static parameters
beta = 1773.6
kappa = 2952.7

#alpha = (.194*(1.7+(1990*aTc**.5588)))/(.2286+(aTc**.5588))
#y =(.0363*4286*(ara**1.775))/(.2286+aTc**1.775)

#initial value
y_scale = .0028
y = convertAraToY(ara[0])
#alpha = convertATCtoAlpha(aTc, y_scale)
alpha = 0

tspan = [0,50]

# ode defininition
def func(t, f):
    x = f
    # differential equations
    dXdt = alpha+beta*(x-1-y+np.sqrt((x+y+1)*(x+y+1)-4*x*y))/(kappa+x-1-y+np.sqrt((x+y+1)*(x+y+1)-4*x*y))-x
    dfdt = [dXdt]
    return dfdt
# test the code above
def tester():
    # solves
    x0_range=[10,100]
    x0= x0_range[0]
    lines=[]
    while x0 < x0_range[1]:
        sol = solve_ivp(func, tspan, [x0])
        ode_line, = plt.plot(sol.t,sol.y.T, label='initial= '+str(x0))
        lines.append(ode_line)
        #lines.append(c_line)
        x0+=20
    plt.legend(handles=lines)
    # display plot
    plt.title('sigma abundance')
    plt.xlabel('Time')
    plt.ylabel('[sigma factor]')
    plt.show()


