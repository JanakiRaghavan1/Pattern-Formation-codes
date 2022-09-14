%Single cell Fitzhugh Nagumo model
clear
%initial condition
v=0.8;%initial fast variable;it can be above or below threshold and weget the action potential accordingly
w=0.;%initial slow variable value
k = 0.6;
a = 0.139;
epsilon = 0.001;
b=0.00;%0.14
t=5e3;
del_t=1;


for i=1:ceil(t/del_t),
%running Fitzhugh Nagumo equation
%if i<ceil(.5/del_t)
   % I=0.01;
%else
    %I=0;
%end
del_v=v.*(1-v).*(v-a)-w;
del_w=epsilon*(k*v-w-b);

%Euler method
v=v+del_v*del_t;
w=w+del_w*del_t;
%increment the values to the new set of values 

vf(i)=v;
wf(i)=w;
end;
 %vf1=vf(1,2001:23580);
%  %wf1=wf(1,2001:23580);
%  x=linspace(0,1,1000);
%  y= x.*(1-x).*(x-a);
%  y1=k*x-b;
%  
%  figure
%  plot(x,y)
% hold on
% plot(x,y1,'r')
% hold on
% plot(vf,wf)
%  hold on
plot(vf)
%plot(vf(1,1253:23580),wf(1253:23580))
%save(sprintf('A_%1.2f_B_%1.2f_singlecell.mat',A,B));