# -----------------------------------------------------------
# demonstrates how to model the original matlab code
# email kdnichols@gmail.com
# This project was created for a class at UW-Madison
# -----------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# parameters
gH = 1.134
gL = 1.188
#r =[.01,1]
r = 1
#f =[.01,1]
f =1

def stochastic_switching(y, t):

    NL,NH = y # unpack y
    dNLdt = -r*NL + f*NH + gL*NL
    dNHdt =  r*NL - f*NH + gH*NH
    return dNLdt,dNHdt # repack dydt


def plot_phase(ODE, u_range, v_range, n_grid=100):
    ax = plt.axes()
    u = np.linspace(u_range[0], u_range[1], n_grid)
    v = np.linspace(v_range[0], v_range[1], n_grid)
    uu, vv = np.meshgrid(u, v)
    u_vel = np.empty_like(uu)
    v_vel = np.empty_like(vv)
    for i in range(uu.shape[0]):
        for j in range(uu.shape[1]):
            u_vel[i,j], v_vel[i,j] = ODE(np.array([uu[i,j], vv[i,j]]), None)

    speed = np.sqrt(u_vel**2 + v_vel**2)
    lw = 0.5 + 2.5 * speed / speed.max()
    ax.streamplot(uu, vv, u_vel, v_vel, linewidth=lw, arrowsize=1.2, density=1)#, color='thistle')
    ax.set_xlabel('N_L')
    ax.set_title('Stochastic Switching phase plot')
    ax.set_ylabel('N_H')
    return ax
def tester():
    #plt.subplot(1, 2, 1)
    ax = plot_phase(stochastic_switching, (0,10), (0,10))

    '''
    NL0 = .2
    NH0 = .2
    plt.subplot(1, 2, 2)
    tspan = np.linspace(0,20)
    sol = odeint(stochastic_switching, tspan, [NL0,NH0])
    print(sol)
    plt.plot(tspan, sol[0]/sol[1])
    plt.xlabel('time')
    plt.ylabel('switching rate')
    '''
    plt.show()

