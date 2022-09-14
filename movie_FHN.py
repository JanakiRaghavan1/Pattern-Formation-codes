from sage.all import *
import os
import sys
import pylab as plt
from FHN import SingleOscillator, Torus
import time
import os.path

assert len( sys.argv ) == 5, 'Usage: sage <program> <b> <D_u> <D_v> <n_x>'

b          = float( sys.argv[1] ) #0.15
D_u        = float( sys.argv[2] ) #0.001
D_v        = float( sys.argv[3] ) #0.001
n_x        = float( sys.argv[4] )
n_y        = n_x
epsilon    = 0.001;
k          = 0.6;
alpha      = 0.139;

dt = 100

#num_transient_frames = 20 # number of transient frames
#num_obs = 20  # number of frames

num_transient_frames = 100000 # number of transient frames
num_obs = 1000  # number of frames

#num_transient_frames = 99000 # number of transient frames
#num_obs = 1000  # number of frames

print('# n_x = %d, n_y = %d, b = %g, D_u = %g, D_v = %g' %(n_x, n_y, b, D_u, D_v))

#-----------Running the system for a single oscillator-----------#
#--------to get a number of index values on the limit cycle------#
T = ode_solver()
T.function = SingleOscillator(epsilon=epsilon, k=k, b=b, alpha=alpha)
T.ode_solve(y_0=[0.1,0.1], t_span=[0,10000], num_points=10) # transience
T.ode_solve(y_0 = T.solution[-1][1], t_span=[0,10000], num_points=10000)
y_single_cell = plt.array( flatten(T.solution)).reshape((10001,3))
#----------------------------------------------------------------#

T.function = Torus(n_x=n_x, n_y=n_y, epsilon=epsilon, k=k, b=b, alpha=alpha, D_u=D_u, D_v=D_v)

num_osc = n_x * n_y
indices = plt.randint(10**4, size=num_osc);
y = list( plt.r_[y_single_cell[indices,1], y_single_cell[indices,2]] ) # Choose random values from LC

T.ode_solve(y_0=y, t_span=[0,dt*num_transient_frames], num_points=10) # Discard transience
T.ode_solve(y_0=T.solution[-1][1], t_span=[0,dt*num_obs], num_points=num_obs) # Data to be saved

i1 = 1
ftest = 0
while ftest==0:
    filename = "FHN_N_%d_b_%06.4f_Du_%06.4f_Dv_%06.4f_%d.mat" %(n_x, b, D_u, D_v, i1)
    if os.path.isfile(filename) == True:
        i1 = i1 + 1
    else:
        ftest = 1

from scipy.io.matlab.mio import savemat
savemat( filename,
        {"n_x" : n_x, "dt" : dt, "num_obs" : num_obs, "num_transient_frames" : num_transient_frames,
            "data" : plt.array( flatten(T.solution) )} )
