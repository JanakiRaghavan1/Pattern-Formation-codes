# Code to run in parallel on annapurna and others
#
# Rajeev Singh
# 20130803: first version
#
# 20130903: second version (Modified by SNM and Janaki)

from sage.all import *
import sys
import pylab as plt

num_osc    = 20
B          = float( sys.argv[1] )
D_u        = float( sys.argv[2] )
D_v        = float( sys.argv[3] )
num_trials = int(   sys.argv[4] )

# We require A in the range (A_hopf, A_hopf+4)
A_in       = B*B + 1
A_fin      = (B*B + 1) + 4
dA         = 0.04

filename = 'B_%1.1f_Dv_%1.2f.dat' %(B, D_v)
f = open(filename, 'w')

print '# num_osc = %g, B = %g, D_u = %g, D_v = %g, num_trials = %g' %(num_osc, B, D_u, D_v, num_trials)

@parallel('multiprocessing', ncpus=8) # number of cpus that will be used
def func(num_osc, A, B, D_u, D_v, num_trials):
    import pylab as plt
    from Brusselator import Brusselator
    a = Brusselator(num_osc=num_osc, A=A, B=B, D_u=D_u, D_v=D_v)
    pattern_count = plt.zeros(8)
    pattern_count[0] = A
    for ii in range(num_trials):
        a.set_initial_condition()
        a.integrate(10000, 10) # transience
        a.integrate(100, 1000)
        label, name = a.which_pattern()

        # If chimera, run the program for a bit longer to confirm
        # (we needn't choose a large num_obs, as we just want to
        #  check if at least one node is not oscillating)
        if label==3:
            a.integrate(1000, 1000)
            label, name = a.which_pattern()

        pattern_count[label+1] += 1 # index 0 is for value of A
    
    return pattern_count

inputs = []
for A in plt.linspace(A_in+dA, A_fin, 100):
    inputs.append( (num_osc, A, B, D_u, D_v, num_trials) )

#print func(num_osc, A, B, D_u, D_v, num_trials)

outputs = func(inputs)

for input, output in outputs:
    sys.stdout = f
    print '% 8.4f' %output[0],
    for ii in range(1,8):
        print '% 4d' %output[ii],
    print
    sys.stdout.flush()
