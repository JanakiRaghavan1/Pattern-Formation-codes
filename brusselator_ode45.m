% A = 3; B = 1; D = 0.1;  %[GS]
% A = 12.0; B = 3.0; D = 0.20; % [GS-APS]
 %A = 6; B = 2; D = 1;    %[SPOD]
 %A = 6; B = 2; D = 0.44; %[CHIMERA]
% A =  7.2; B = 2.2; D = 0.17; % fig2 [Single defect; also APS in GS]
%A =  4.3; B = 0.7; D = 0.41; % [BURSTS]
% A =  9.5; B = 2.2; D = 0.22; % Here we get travelling chimera
% A =  5.3; B = 2.0; D = 0.20; % newfig9 [APS in GS (second)]
% A =  5.0; B = 0.8; D = 0.25; % [GS + defect (mismatch)]
% A = 10.0; B = 2.9; D = 0.13; % newfig2 [APS + defect]
% A =  7.0; B = 2.1; D = 0.28; % region between GS and SPOD [SPOD + defect]
%A=5.5;B=2,D=0.2;
%A =  6.0; B = 2.0; D = 0.20; % newfig7 [APS + defect (other); also APS]
 %A =  9.5; B = 2.2; D = 0.236; % Here we get travelling SPOD
%%%% A =  7.7; B = 2.3; D = 0.25; % newfig 3
 %A =  7.2; B = 2.3; D = 0.30; % newfig4
%A = 10.2; B = 2.3; D = 0.20; % newfig5
% A =  7.3; B = 2.2; D = 0.26; % newfig6
% A =  7.2; B = 2.2; D = 0.19; % fig2 old not frequent
%%%% A =  5.8; B = 2.0; D = 0.20; % newfig8: also APS with defect
 B = 0.6; A=2+B*B; D = 0.8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Paperdata
%A=4;B=0.9;D=0.25;%Figure1)a1 %defect in GS
%A=9.5;B=2.2;D=0.237;%Figure1)a2 %travelling chimera
%A=4;B=1;D=0.8;%Figure1)a3 %shuffling
%A=6;B=2.;D=0.2; %APS
%A=5.3;B=2;D=0.2;%Figure1)a4%APS and GS combination
%A=10;B=2.9;D=0.13;%Figure1)a5^APS defect
%A=7;B=2.1;D=0.28;%spod+defect
D = D/2;
D0 = D;

N = 40; % number of oscillators

U = rand(1,2*N);

opt  = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
opt2 = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

%============Obtain values from limit cycle=============%
params = [A B 0];                                       %
                                                        %
T0 = 0:1e2:1e3;                                         %
T1 = 0:0.005:10;                                        %
[t, UL0] = ode45(@brusselator_sys, T0, U, opt, params); %
[t, UL1] = ode45(@brusselator_sys, T1, UL0(end,:), opt, params);
                                                        %
u0 = UL1(:,1:N);                                        %
v0 = UL1(:,N+1:2*N);                                    %
                                                        %
fprintf('Limit cycle for D=0 has been obtained\n');     %
%=======================================================%

params = [A B D];
T = 0:0.01:350;

for ii = 1 : 1
    
    fprintf('(Run %02d) Obtaining data\n',ii);
    
    for jj = 1:N
        rndindx = ceil(1e3*rand(1));%randi(1000);
        u_in(jj) = u0( rndindx, 1);
        v_in(jj) = v0( rndindx, 1);
    end
    
    U0 = [u_in v_in];
    
    [t, U1] = ode45(@brusselator_sys, T0, U0, opt, params);
    [t, UV] = ode45(@brusselator_sys, T, U1(end,:), opt2, params);

    u = UV(:, 1:N);
    v = UV(:, N+1:2*N);
    
    
    %clf;
    figure, imagesc(1:40, t, u);
     figure,plot(t, u); %pause(1e-1);
    
    %save(sprintf('A_%1.2f_B_%1.2f_D_%1.2f_%d.mat',A,B,D,ii),'u','v','t');
    
end

clear T U UV params opt u0 v0 u_in v_in r1 r2 ii jj;