from sage.all import *
import os
import sys
import pylab as plt
from Brusselator import Brusselator
import time
import os.path

assert len( sys.argv ) == 6, 'Usage: sage <program> <A> <B> <D_u> <D_v> <n_x>'

A          = float( sys.argv[1] )
B          = float( sys.argv[2] )
D_u        = float( sys.argv[3] )
D_v        = float( sys.argv[4] )
n_x        = float( sys.argv[5] )
n_y        = n_x

print '# num_osc = %g, A = %g, B = %g, D_u = %g, D_v = %g' %(n_x, A, B, D_u, D_v)

#@parallel('multiprocessing', ncpus=8) # number of cpus that will be used

a = Brusselator(num_osc=n_x*n_y, n_x=n_x, A=A, B=B, D_u=D_u, D_v=D_v, bc='Torus')

a.set_initial_condition()
a.integrate(100000, 10) # transience (normal: 100000)
a.integrate(100, 5000)

i1 = 1
ftest = 0
while ftest==0:
    filename = "BR_N_%d_A_%06.4f_B_%06.4f_Du_%06.4f_Dv_%06.4f_%d.mat" %(n_x, A, B, D_u, D_v, i1)
    if os.path.isfile(filename) == True:
        i1 = i1 + 1
    else:
        ftest = 1

from scipy.io.matlab.mio import savemat
savemat( filename,
        {"n_x" : n_x, "dt" : 0, "num_obs" : 0, "num_transient_frames" : 0,
            "data" : plt.array( flatten(a.T.solution) )} )
