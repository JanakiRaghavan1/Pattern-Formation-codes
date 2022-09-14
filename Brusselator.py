# Python class to simulate brusselator periodic ring using Sage's
# wrapper to GSL.
#
# Rajeev Singh
# 20130803: first version
#
# 20130903: second version (Modified by SNM and Janaki)

from sage.all import *
from brusselator_rhs import PeriodicRing
import pylab as plt
import numpy as np
from numpy import *
from scipy.signal import argrelextrema

class Brusselator:
    def __init__(self, num_osc, A, B, D_u, D_v):
        self.num_osc = num_osc
        self.A       = A
        self.B       = B
        self.D_u     = D_u
        self.D_v     = D_v

        self.T = ode_solver()
        self.T.algorithm = 'rkf45' # other GSL algorithms can be used here
        # setting the integrator to get the unperturbed limit cycle
        # for given parameters
        self.T.function = PeriodicRing(num_osc=2, A=A, B=B, D_u=0, D_v=0)
        self.T.ode_solve(y_0=[0.1,0.1,0.1,0.1], t_span=[0,100000], num_points=10) # transience
        self.T.ode_solve(y_0 = self.T.solution[-1][1], t_span=[0,1000], num_points=10000)
        self.y_single_cell = plt.array( flatten(self.T.solution)).reshape((10001,5))

        self.T.function = PeriodicRing(num_osc=num_osc, A=A, B=B,
                D_u=D_u, D_v=D_v) # setting up the given system
    
        # setting up random initial condition
        if plt.std( self.y_single_cell[:,1] ) < 1e-3 :
            self.isOscillating = False
        else :
            self.isOscillating = True

        self.initial_condition = None
        self.set_initial_condition()

        # variable to hold the trajectory
        self.trajectory = None

    def set_initial_condition(self):
        if self.isOscillating:
            num_obs_single_cell = self.y_single_cell.shape[0]
            indices = plt.randint(num_obs_single_cell, size=self.num_osc)
            # since to get unperturbed solution we have simulated two
            # uncoupled oscillators, we use the indices (1,3) for the
            # (u,v) component of the solution for the first oscillator.
            self.initial_condition = list( plt.r_[self.y_single_cell[indices,1], self.y_single_cell[indices,3]] )
        else : # return random numbers if the system is not oscillating
            self.initial_condition = list( plt.random(2*self.num_osc) )

    def integrate(self, t_total=100, num_obs=10):
        self.T.ode_solve(y_0=self.initial_condition,
                t_span=[0,t_total], num_points=num_obs-1)
        # resetting initial condition to the updated state
        self.initial_condition = self.T.solution[-1][1]
        # saving the trajectory
        self.trajectory = plt.array( flatten(self.T.solution)).reshape((num_obs,2*self.num_osc+1))
    
    def which_pattern(self, is_oscillating_cutoff=1e-4, is_sync_cutoff=1e-5):
        num_osc = self.num_osc
        y = self.trajectory
        num_obs = y.shape[0]

        # checking for AD
        if plt.var( y[:,1:num_osc+1] ) < is_oscillating_cutoff:
            return 5, 'AD'

        # calculating the amplitude of each oscillator
        amplitude = plt.var( y[:,1:num_osc+1], axis=0 )

        # counting the number of non-oscillating elements
        num_non_oscillating = plt.sum( amplitude < is_oscillating_cutoff )

        # checking for SPOD and chimera
        if num_non_oscillating == num_osc:
            return 2, 'SPOD'
        elif num_non_oscillating > 0 :
            return 3, 'chimera'

        # checking for SO
        if plt.mean( plt.var( y[:,1:num_osc+1],axis=1 ) ) < is_sync_cutoff:
            return 0, 'SO'

        # checking for APS by lookng whether odd oscillators are inphase
        if plt.mean( plt.var( y[:,1:num_osc+1:2  ],axis=1 ) ) < is_sync_cutoff:
            return 1, 'APS'

        #------------------------------------------------------------------------------#
        #                      Check for gradient synchronization                      #
        #------------------------------------------------------------------------------#
        # Initialization
        GS_var  = 0
        GS_diff = 0
        osc_rang = np.arange(0,num_osc)

        GS_diff_test = list(osc_rang) # Check GS for each consecutive oscillator pair
        uindx   = list(osc_rang)      # Indices of maxima of all oscillators
        ulen    = list(osc_rang)      # Number of maxima seen for each oscillator
        ypk_var = list(osc_rang)      # Variance of the values at each peak

        for i1 in osc_rang:
            get_peakindx = argrelextrema(y[:,i1+1], np.greater)
            uindx[i1] = get_peakindx[0]
            ulen[i1]  = len(uindx[i1])

        Lm = min(ulen)

        # If there are less than three peaks, then this is not GS
        if Lm < 3:
            return 6, 'others'

        for j1 in osc_rang:
            u0 = uindx[j1][1:Lm-1]
            u1 = uindx[(j1 + 1) % 20][1:Lm-1]
            ypk_var[j1] = plt.var(y[ u0[1:Lm-1], 1 ])

            udiff = u0 - u1
            if len(udiff) == 0:
                udiff = plt.zeros(2)

            # If GS, then the signs of udiff are either all +ve or all -ve
            GS_diff_test[j1] = 0
            if max(udiff) <= 0 or min(udiff) >= 0:
                GS_diff_test[j1] = 1

        # Is the GS test true for all oscillator pairs?
        if min(GS_diff_test)==1:
            GS_diff = 1

        # If the variance of the values at the peaks is large, it is not GS
        uvar = max(ypk_var)
        if uvar<5e-2:
            GS_var = 1        

        #------------------------------------------------------------------------------#

        # If the above is satisfied, we register the behaviour as GS
        if GS_diff == 1 and GS_var == 1:
            return 4, 'GS'

        # if no pattern among the above is detected return 'others'
        return 6, 'others'

    def plot(self):
        y = self.trajectory
        plt.figure()
        plt.pcolor( plt.r_[0:y.shape[1]/2+1], y[:,0], y[:,1:self.num_osc+1] )
        plt.colorbar()
        plt.xlabel('oscillator index')
        plt.ylabel('time')
        plt.axis('tight')
        plt.savefig('save.png')

if __name__ == '__main__':
    #a = Brusselator(num_osc=20, A=5.5, B=2, D_u=0.0, D_v=0.2) # Here we get APS or others
    #a = Brusselator(num_osc=20, A=3, B=1, D_u=0.0, D_v=0.05)  # Here we get SO or GS
    #a = Brusselator(num_osc=20, A=4.5, B=1, D_u=0.0, D_v=0.2)  # Here we get SO or GS
    a = Brusselator(num_osc=20, A=5.2, B=2, D_u=0.0, D_v=0.2)  # Here we get SO or GS
    
    indx=0
    while (indx!=4 and indx!=6):
        a.set_initial_condition()
        a.integrate(10000, 10) # transience
        a.integrate(200, 2000)
        indx, p = a.which_pattern()
        print p
        # If chimera, run the program for a bit longer to confirm
        # (we needn't choose a large num_obs, as we just want to
        #  check if at least one node is not oscillating)
        if indx==3:
            a.integrate(1000, 1000)
            indx, p = a.which_pattern()
            print p
    a.plot()
