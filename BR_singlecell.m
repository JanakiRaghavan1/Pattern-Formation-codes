%Single cell Fitzhugh Nagumo model
clear
%initial condition
v=0.0;%initial fast variable;it can be above or below threshold and weget the action potential accordingly
w=0;%initial slow variable value
A=4.3;
B=0.7;
t=5e2;
del_t=0.01;


for i=1:ceil(t/del_t),
%running Fitzhugh Nagumo equation
%if i<ceil(.5/del_t)
   % I=0.01;
%else
    %I=0;
%end
del_v=B - (A+1)*v + (v.^2).*w;
del_w=(A*v - v.^2.*w);

%Euler method
v=v+del_v*del_t;
w=w+del_w*del_t;
%increment the values to the new set of values 

vf(i)=v;
wf(i)=w;
end;
 %vf1=vf(1,2001:23580);
 %wf1=wf(1,2001:23580);
 plot(vf)
%plot(vf(1,1253:23580),wf(1253:23580))
%save(sprintf('A_%1.2f_B_%1.2f_singlecell.mat',A,B));