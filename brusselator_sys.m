function dydt = brusselator_sys(t,U,params)

A = params(1);
B = params(2);
D = params(3);

N = size(U,1)/2;
dydt = zeros(N,1);

for n1 = 1:N
    u = U(0*N + n1);
    v = U(1*N + n1);
    
    if (n1>1)
        vm = U( 1*N  + n1 - 1);
    else
        vm = U( 2*N );
    end
    if (n1<N)
        vp = U( 1*N + n1 + 1);
    else
        vp = U( 1*N +1);
    end
    
    dydt(0*N + n1) = B - (A+1)*u + (u^2)*v;
    dydt(1*N + n1) = (A*u - u^2*v) + D*(vp + vm - 2*v);
end

end