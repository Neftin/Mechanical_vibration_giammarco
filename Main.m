
clear all

%try to impement the lsqnonlin within the simulation for the output!!!

load('dataMass');

rng default 

rg = (1:length(t));

global v t x1 x2 x3 gain_x;

gain_x = (0.0706/16000); %encoder!!!

%data processing and homologation, and choice of range in samples

INPUT_FOR_ESTIMATION = timeseries(v(rg),t(rg));

OUTPUT_TO_TEND       = timeseries([x1(rg)*gain_x , x2(rg)*gain_x, x3(rg)*gain_x],t(rg));

comparison_data      = iddata([x1(rg)*gain_x , x2(rg)*gain_x , x3(rg)*gain_x], v(rg),Ts);

%initial guess

%K = [k1 k2 k3];

%P = [ A_mean B1(3:7) B2(4:7) B3(5:7) ]; %first guess by step response estimation (better if manually)

%the zeroes must be dropped!

%first guess

P = [1.5 1.5 1.5 1 1 1 1 1 5];





%P = [m1,m2,m3,c1,c2,c3,c12,c23,g_v]; %MANAGE LOWERBOUNDs AND UPPERBOUNDS, k
%given
% 
 lb = zeros(1,9).';
 ub = [2.1 2.1 2.1 10 10 10 10 10 7];

%% lsq non lin full damping

options = optimoptions('lsqnonlin','Display','iter','StepTolerance',0.000000001,'OptimalityTolerance',0.0000001);

param = lsqnonlin(@error_statespace_full,P,lb,ub,options)

%% sys reconstruction

m1      = param(1);
m2      = param(2);
m3      = param(3);

c12     = param(4);
c23     = param(5);

c1      = param(6);
c2      = param(7);
c3      = param(8);

g_v_est_full = param(9);

I = eye(3);

Z = zeros(3);

M = [m1 0 0;  
     0 m2 0; 
     0 0 m3];
 
 %K definition fixed
 
 k1=774;k2=770;k3=396;
 
 K = [k1   -k1        0;  
     -k1 k1+k2      -k2; 
      0    -k2     k2+k3];
  
 Cc = [+c1+c12   -c12        0
      -c12  +c2+c12+c23     -c23
        0        -c23      +c3+c23] %called C but it isn't a problem in order of assignments
 
 b = [g_v_est_full 0 0].';
 
 A = [Z I; -M\K -M\Cc];    %left divide for the inverse
 B = [Z(:,1); M\b]; %single input
 C = [I Z];               %multiple output
 D = Z(:,1);

%%

sys_guess_full = ss(A,B,C,D);

figure(1);
compare(comparison_data,sys_guess_full)

%% lsq non lin proportional damping

P = [1.5 1.5 1.5 1.5 0.0001 6 ];
 lb = zeros(1,9).';
 ub = [2.1 2.1 2.1 10 10 10 10 10 7];

options = optimoptions('lsqnonlin','Display','iter','StepTolerance',0.000000001,'OptimalityTolerance',0.0000001);

param = lsqnonlin(@error_statespace_proportional,P,lb,ub,options)

%% sys reconstruction

m1      = P(1);
m2      = P(2);
m3      = P(3);

ca      = P(4);
cb      = P(5);

g_v_est_prop = P(6);





I = eye(3);

Z = zeros(3);

M = [m1 0 0;  
     0 m2 0; 
     0 0 m3];
 
 %K definition fixed
 
 k1=774;k2=770;k3=396;
 
 K = [k1   -k1        0;  
     -k1 k1+k2      -k2; 
      0    -k2     k2+k3];
  
 C = ca*M + cb*K; %called C but it isn't a problem in order of assignments
    
 
 
 b = [g_v_est_prop 0 0].';
 
 A = [Z I; -M\K -M\C];    %left divide for the inverse
 B = [Z(:,1); M\b]; %single input
 C = [I Z];               %multiple output
 D = Z(:,1);
%% comparison

sys_guess_prop = ss(A,B,C,D);

figure(2);
compare(comparison_data,sys_guess_prop)

pp_impulse

%% MODAL ANALYSIS - eigenvalues

%find eigenvalues as diagonal matrix and eigenvectors

[ U_eig omega2_Eig_matrix ] = eig(K,M);

    for i=1:3
        U_Eig(1:3,i) = U_eig(1:3,i)./U_eig(1,i);
        omega_eig(i) = omega2_Eig_matrix(i,i)^(1/2);
    end

clear omega2_Eig_matrix

%% MODAL ANALYSIS - RAILEIGHT

%in this case we can use the data from the proportional since we are
%interested only in the undamped version!

%variation mode shape vector for the optimization (note that the real space
%is an R^2)
%space, so we can plot the surface if we want.

%Raileight quotient
R_q = @(u) ( [1; u(1); u(2)].' * K * [1; u(1); u(2)] )  /  ( [1; u(1); u(2)].' * M * [1; u(1); u(2)] ) ;

R_q([1;1])

%sintax u_R[i] are the real modes, uu_R are only temporary variables
%multiplied times the kernels to solve Raileight minimas


u0_R = [1,1]
%find minima
[uu_R1 , omega1_R ] = fminunc(R_q,u0_R)

%per trovare gli altri modi io introdurrei semplicemente il vincolo di
%ortogonalità in questo modo è una figata perchè vai a tagliare la
%superficie e trovare un minimo su una dimensione in meno!

u_R1 = [1 uu_R1].'; %first mode

B_Ker_u1 = null([1 uu_R1]);

R_q2 = @(uu) ( ( ( [1; uu].'*B_Ker_u1.' ) * K * ( [1; uu].'*B_Ker_u1.' ).' )...
         /     ( ( [1; uu].'*B_Ker_u1.' ) * M * ( [1; uu].'*B_Ker_u1.' ).' ) ) ;

% [ 1 uu has to be multiplied times B_Ker_u
uu0_R = 1;
     
[uu_R , omega2_R ] = fminunc(R_q2,uu0_R)

u_R2 = ([1; uu_R].'*B_Ker_u1.').' 
u_R2 = u_R2.*(1/u_R2(1));    % second mode normalize at the first entry

u_R3 = null([u_R2.' ; u_R1.'])
u_R3 = u_R3.*(1/u_R3(1));

U_Rail     = [u_R1 , u_R2 , u_R3]
omega_Rail = [R_q(u_R1(2:3)) R_q(u_R2(2:3)) R_q(u_R3(2:3))].^(1/2)

pp_raileight

%% MODAL ANALYSIS - Matrix iteration method

%Essentially it is a way to find the largest eigen value and the
%corresponding eigenvector. so the procedur eis to find one by one alla the
%omegas. combining MIM-deflation-MIM-deflation-MIM

    %deflation: now we will move this eigenvalue to zero (theoretically), also
    %if it is not zero no problem because this operatin is for remove it from
    %the pole position of the eigenvalues

    %initialization and frst guess
uMIM      = [1 1 1].';
D_0       = K\M;
U_MIM     = zeros(3);

eig_D_MIM = zeros(3,1);

D = D_0;
for i=1:3
%Iterate for the 3 modes
    for j=1:15
        if (j==1)
            uMIM = ones(3,1); %first guess again
        end
        uMIM = D*uMIM;
        uMIM = uMIM./uMIM(1);
    end
    %compute i-th eigenvalue (SEE PAGE STANDARD EIGENVALUE PROBLEM ON YOUR NOTEBOOK TO PUT COMPUTATIONS ON REPORT)
    eig_D_MIM(i) = (uMIM.' * D * uMIM) / (uMIM.'*uMIM);
    U_MIM(:,i)   = uMIM; 
    %Matrix deflation:
        normalizator = uMIM.' * M * uMIM;     % compute the normalization factor
        uMIM = uMIM./( (normalizator)^(1/2) ); % use the normalizaed-to-mass uMIM vector to deflate the matrix D
    D = D - eig_D_MIM(i)*uMIM*uMIM.'*M;
end

[U_MIM,U_Rail,U_Eig]

%% MODAL ANALYSIS - Laplace

%tranfer function with force (voltage-to-force coefficient compensation)
LaPlace_full = (1/g_v_est_full)*tf(sys_guess_full);

figure(3)
bode(LaPlace_full(1))
grid minor;

%tranfer function with force (voltage-to-force coefficient compensation)
LaPlace_prop = (1/g_v_est_prop)*tf(sys_guess_prop);

figure(4)
bode(LaPlace_prop(1))
grid minor;

%aggiungi post_proc


