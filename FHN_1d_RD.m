% 1D_RD_FHN
% This program simulates the 1D FHN RD equation
% with periodic boundary conditions using
% a Crank-Nicolson method with an LU decomposition (full pivoting)
% and 2-step Adams-Bashforth for the nonlinear terms

clear; clc; clf; %close all;
nt = 5e4;     % number of timesteps
dt = 1;     % timestep
trans = 3*nt/4;

N  = 20;   % grid length
dx = 1; % grid spacing

Dv = 0;
Dw = 0.6e-3;%0.5e-3--APS ,0.6e-3-->GAPS,GS, 1e-3--others,1.5e-3/1.6e-3--SO, 2.5e-3--chimera, 4e-3--SPOD
k = 0.6;
a = 0.139;
epsilon = 0.001;
b= 0.06;%0.06--APS,others , 0.14--GS ,0.07--SO, 0.08--chimera,SPOD
sv = (Dv/2)*dt/(2*(dx^2));
sw = (Dw/2)*dt/(2*(dx^2));
%fprintf('s  = %2.16g',s);  disp(' ');

% Initial condition
V0 = zeros(N,1);
W0 = zeros(N,1);
v = rand; w = rand;
for ii = 1:1e4
     del_v = v.*(1-v).*(v-a)-w;
     del_w = epsilon*(k*v-w-b);
     vf(ii) = v+del_v*dt;
     wf(ii)= w+del_w*dt;
     
      v = vf(ii);
      w = wf(ii);
end
for jj=1:N
    indx = floor(rand*1e4);
    V0(jj) = vf(indx);
    W0(jj) = wf(indx);
end
V1 = V0;W1 = W0;

%  [~, y_single_cell ] = ode15s(@single_cell, 0:500:1.0e4, rand(2,1) );
%  [t, y_single_cell ] = ode15s(@single_cell, 0:1:1.0e4, y_single_cell(end,:) );
%  indices = randi(10^4, N,1);
%  V0 = y_single_cell(indices,1);
%  W0 = y_single_cell(indices,2);
% V0form = 0.1*exp(-20*x.^2);
% V0 = V0form';

o1 = ones(N,1); o2 = ones(N-1,1);
S2v = diag((1+2*sv)*o1) + diag(-sv*o2,1) + diag(-sv*o2,-1);
S1v = diag((1-2*sv)*o1) + diag( sv*o2,1) + diag( sv*o2,-1);
S2w = diag((1+2*sw)*o1) + diag(-sw*o2,1) + diag(-sw*o2,-1);
S1w = diag((1-2*sw)*o1) + diag( sw*o2,1) + diag( sw*o2,-1);

%I = diag(ones(N,1));
%[sum(sum(S1v-I)) sum(sum(S2v-I)) sum(sum(S1w-I)) sum(sum(S2w-I))]
%return

% Periodic BCs:
S2v(1,end) = S2v(1,2); S2v(end,1) = S2v(end,end-1);
S1v(1,end) = S1v(1,2); S1v(end,1) = S1v(end,end-1);
S2w(1,end) = S2w(1,2); S2w(end,1) = S2w(end,end-1);
S1w(1,end) = S1w(1,2); S1w(end,1) = S1w(end,end-1);

SS1v = sparse(S1v);
SS2v = sparse(S2v);
SS1w = sparse(S1w);
SS2w = sparse(S2w);

%LUPQ Decomposition to get inverse
[L_v,U_v,P_v,Q_v] = lu(SS2v);
[L_w,U_w,P_w,Q_w] = lu(SS2w);

% Vmax = zeros(nt,1);
% Vchg = zeros(nt,1);

V_rec = zeros(nt-trans,N);

% Plot the solution profile every pl steps
pl = 1;
for ii=1:nt

     reac_v0 = V0.*(1-V0).*(V0-a)-W0;
     reac_w0 = epsilon*(k*V0-W0-b);
     
     reac_v1 = V1.*(1-V1).*(V1-a)-W1;
     reac_w1 = epsilon*(k*V1-W1-b);
     
     % Two-step Adams-Bashforth
     b_v = SS1v*V1 + dt*((3/2)*reac_v1 - (1/2)*reac_v0);
     b_w = SS1w*W1 + dt*((3/2)*reac_w1 - (1/2)*reac_w0);
     
     V2 = Q_v*(U_v\(L_v\(P_v*b_v)));
     W2 = Q_w*(U_w\(L_w\(P_w*b_w)));
    
     V0 = V1;W0 = W1;
     V1 = V2;W1 = W2;
     
     if ii>trans
     V_rec(ii-trans,:) = V2;
     end

%      if mod(ii,pl)==0
%      imagesc(V2);
%      %title(strcat('1D RD equation: \partial V/\partial t',...
%      %' = \D partial^2 V/\partial x^2 + V(1-V)'));
%      %xlabel('x','FontSize',14);
%      %ylabel('V(x,t)','FontSize',14);
%      %axis([x(1) x(end) -0.1 1.1]);
%      %text(-9,-0.05,sprintf('chg = %2.16f',Vchg(ii)));
%      pause(0.1);
%      end
end
imagesc(V_rec);
set(gca,'fontsize',20)
 xlabel('Oscillator #','FontSize',14);
ylabel('V(x,t)','FontSize',14);
%plot(V_rec)
