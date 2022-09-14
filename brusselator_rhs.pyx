


# 
# Cython code to calculate the RHS for brusselator ring so that the
# GSL wrapper of Sage can be used for simulation.
#
# NOTE: We explicitly divide the diffusion terms by 2 (number of
# neighbors) so that the results can be directly compared to that with
# arbitrary coupling topology.
#
# Rajeev Singh
# 20130803: first version
#

cimport sage.gsl.ode
import  sage.gsl.ode
include 'gsl.pxi'

cdef class PeriodicRing(sage.gsl.ode.ode_system):
    cdef int    num_osc
    cdef double A, B
    cdef double D_u, D_v

    def __cinit__(self, int num_osc, double A, double B, double D_u, double D_v):
        self.num_osc = num_osc
        self.A       = A
        self.B       = B
        self.D_u     = D_u
        self.D_v     = D_v
        
    cdef int c_f(self, double t, double *y, double *dydt):
        cdef int i
        cdef int n = self.num_osc
        # pointers to the input and output arrays for convenience
        cdef double *u    = &y[0]
        cdef double *v    = &y[n]
        cdef double *dudt = &dydt[0]
        cdef double *dvdt = &dydt[n]

        # left site
        dudt[0] = self.B + u[0]*u[0]*v[0] - u[0] - self.A*u[0] \
                  + self.D_u * (u[n-1] + u[1] - 2*u[0] ) / 2
        dvdt[0] = self.A*u[0] - u[0]*u[0]*v[0]  \
                  + self.D_v * (v[n-1] + v[1] - 2*v[0] ) / 2
                
        # bulk
        for i in range(1,n-1):
            dudt[i] = self.B + u[i]*u[i]*v[i] - u[i] - self.A*u[i] \
                      + self.D_u * (u[i-1] + u[i+1] - 2*u[i] ) / 2
            dvdt[i] = self.A*u[i] - u[i]*u[i]*v[i]  \
                      + self.D_v * (v[i-1] + v[i+1] - 2*v[i] ) / 2
        # right site
        dudt[n-1] = self.B + u[n-1]*u[n-1]*v[n-1] - u[n-1] - self.A*u[n-1] \
                    + self.D_u * (u[n-2] + u[0] - 2*u[n-1] ) / 2
        dvdt[n-1] = self.A*u[n-1] - u[n-1]*u[n-1]*v[n-1]  \
                    + self.D_v * (v[n-2] + v[0] - 2*v[n-1] ) / 2
            
        return GSL_SUCCESS
