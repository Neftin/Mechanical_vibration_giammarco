
clear all
close all

%% SS ANAL


step_response_analysis

%% start
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

P = [1.5 1.4 1.2 0.01 0.001 3 2 2 6.3];

%P = [m1,m2,m3,c1,c2,c3,c12,c23,g_v]; %MANAGE LOWERBOUNDs AND UPPERBOUNDS, k
%given
% 
 lb = zeros(1,9).';
 ub = [2.1 2.1 2.1 10 10 10 10 10 7];
 
%% lsq non lin full damping

options = optimoptions('lsqnonlin','Display','iter','StepTolerance',0.00000000000001,'OptimalityTolerance',0.0000000000001);

[ param_f, ~, residual_f ] = lsqnonlin(@error_statespace_full,P,lb,ub,options);

%% sys reconstruction

m1_f      = param_f(1);
m2_f      = param_f(2);
m3_f      = param_f(3);

c12_f     = param_f(4);
c23_f     = param_f(5);

c1_f      = param_f(6);
c2_f      = param_f(7);
c3_f      = param_f(8);

g_v_est_f = param_f(9);

I = eye(3);

Z = zeros(3);

M_f = [m1_f 0 0;  
     0 m2_f 0; 
     0 0 m3_f];
 
 %K definition fixed
 
 k1=774;k2=770;k3=396;
 
 K_f = [k1   -k1        0;  
     -k1 k1+k2      -k2; 
      0    -k2     k2+k3];
  
 Cc_f = [+c1_f+c12_f   -c12_f        0
      -c12_f  +c2_f+c12_f+c23_f     -c23_f
        0        -c23_f      +c3_f+c23_f]; %called C but it isn't a problem in order of assignments
 
 b_f = [g_v_est_f 0 0].';
 
 A_f = [Z I; -M_f\K_f -M_f\Cc_f];    %left divide for the inverse
 B_f = [Z(:,1); M_f\b_f]; %single input
 C_f = [I Z];             %multiple output
 D_f = Z(:,1);

%% COMPARISON AND NRMSE

sys_guess_full = ss(A_f,B_f,C_f,D_f);



%% lsq non lin proportional damping

P = [1.5 1.5 1.5 1.5 0.0001 6 ];
lb = zeros(1,9).';
ub = [2.1 2.1 2.1 10 10 10 10 10 7];

options = optimoptions('lsqnonlin','Display','iter','StepTolerance',0.000000001,'OptimalityTolerance',0.0000001);

param_p = lsqnonlin(@error_statespace_proportional,P,lb,ub,options)

%% sys reconstruction

m1_p      = param_p(1);
m2_p      = param_p(2);
m3_p      = param_p(3);

ca_p      = param_p(4);
cb_p      = param_p(5);

g_v_est_p = param_p(6);

I = eye(3);

Z = zeros(3);
1
M_p = [m1_p 0 0;  
       0 m2_p 0; 
       0 0 m3_p];
 
 %K definition fixed
 
 k1=774;k2=770;k3=396;
 
 print2file(k1,'report\result\','%3.0f','\n')
 print2file(k2,'report\result\','%3.0f','\n')
 print2file(k3,'report\result\','%3.0f','\n')
 
 K_p = [k1   -k1        0;  
     -k1 k1+k2      -k2; 
      0    -k2     k2+k3];
  
 Cc_p = ca_p*M_p + cb_p*K_p; %called C but it isn't a problem in order of assignments
    
 
 
 b_p = [g_v_est_p 0 0].';
 
 A_p = [Z I; -M_p\K_p -M_p\Cc_p];    %left divide for the inverse
 B_p = [Z(:,1); M_p\b_p]; %single input
 C_p = [I Z];               %multiple output
 D_p = Z(:,1);
%% comparison

sys_guess_prop = ss(A_p,B_p,C_p,D_p);

figure(2);
compare(comparison_data,sys_guess_prop);


pp_estimation
%% MODAL ANALYSIS - eigenvalues - full

%find eigenvalues as diagonal matrix and eigenvectors

[ U_eig_f omega2_eig_matrix_f ] = eig(K_f,M_f);

    for i=1:3
        U_eig_f(1:3,i) = U_eig_f(1:3,i)./U_eig_f(1,i);
        w_f(i) = omega2_eig_matrix_f(i,i)^(1/2);
    end

clear omega2_eig_matrix_f

%% MODAL ANALYSIS - eigenvalues - prop

%find eigenvalues as diagonal matrix and eigenvectors

[ U_eig_p omega2_eig_matrix_p ] = eig(K_p,M_p);

    for i=1:3
        U_eig_p(1:3,i) = U_eig_p(1:3,i)./U_eig_p(1,i);
        w_p(i) = omega2_eig_matrix_p(i,i)^(1/2);
    end

clear omega2_eig_matrix_p

pp_eig

%% MODAL ANALYSIS - RAILEIGHT - full il fallimento [NON FUNZIONA]
% 
% %variation mode shape vector for the optimization (note that the real space
% %is an R^2)
% %space, so we can plot the surface if we want.
% 
% %Raileight quotient
% R_q_f = @(u) ( [1; u(1); u(2)].' * K_f * [1; u(1); u(2)] )  /  ( [1; u(1); u(2)].' * M_f * [1; u(1); u(2)] ) ;
% 
% %sintax u_R[i] are the real modes, uu_R are only temporary variables
% %multiplied times the kernels to solve Raileight minimas
% 
% u0_R = [1,1]; %first guess
% %find minima
% [uu_R1 , omega1_R ] = fminunc(R_q_f,u0_R); %%%%
% 
% %per trovare gli altri modi io introdurrei semplicemente il vincolo di
% %ortogonalità in questo modo è una figata perchè vai a tagliare la
% %superficie e trovare un minimo su una dimensione in meno!
% 
% u_R1 = [1 uu_R1].'; %first mode %%%%
% 
% B_Ker_u1 = null(U_eig_f(:,1).'); %%%%
% 
% R_q2_f = @(uu) ( [uu(1) uu(2)]*B_Ker_u1.' * K_f * ([uu(1) uu(2)]*B_Ker_u1.').' ) /...
%                ( [uu(1) uu(2)]*B_Ker_u1.' * M_f * ([uu(1) uu(2)]*B_Ker_u1.').' );
% 
% % [ 1 uu ] has to be multiplied times B_Ker_u
% uu0_R = [1 1];
%      
% [uu_R2 , omega2_R ] = fminunc(R_q2_f,uu0_R);
% 
%  u_R2 = ([uu_R2]*B_Ker_u1.').' 
%  u_R2 = u_R2.*(1/u_R2(1));    % second mode normalize at the first entry
% 
% u_R3 = null([u_R2.' ; u_R1.']);
% u_R3 = u_R3.*(1/u_R3(1));
% 
% U_Rail_f     = [u_R1 , u_R2 , u_R3];
% omega_Rail_f = [R_q_f(u_R1(2:3)) R_q_f(u_R2(2:3)) R_q_f(u_R3(2:3))].^(1/2);
% 
% %pp_raileight
% 
% [U_Rail_f U_eig_f]
% 
% %% MODAL ANALYSIS - RAILEIGHT - prop [NON FUNZIONA]
% 
% %variation mode shape vector for the optimization (note that the real space
% %is an R^2)
% %space, so we can plot the surface if we want.
% 
% %Raileight quotient
% R_q_p = @(u) ( [1; u(1); u(2)].' * K_p * [1; u(1); u(2)] )  /  ( [1; u(1); u(2)].' * M_p * [1; u(1); u(2)] ) ;
% 
% %sintax u_R[i] are the real modes, uu_R are only temporary variables
% %multiplied times the kernels to solve Raileight minimas
% 
% u0_R = [1,1]; %first guess
% %find minima
% [uu_R1 , omega1_R ] = fminunc(R_q_p,u0_R);
% 
% %per trovare gli altri modi io introdurrei semplicemente il vincolo di
% %ortogonalità in questo modo è una figata perchè vai a tagliare la
% %superficie e trovare un minimo su una dimensione in meno!
% 
% u_R1 = [1 uu_R1].'; %first mode
% 
% B_Ker_u1 = null(u_R1.');
% 
% R_q2_p = @(uu) ( ( ( [1; uu].'*B_Ker_u1.' ) * K_p * ( [1; uu].'*B_Ker_u1.' ).' )...
%          /       ( ( [1; uu].'*B_Ker_u1.' ) * M_p * ( [1; uu].'*B_Ker_u1.' ).' ) ) ;
% 
% % [ 1 uu ] has to be multiplied times B_Ker_u
% uu0_R = 1;
%      
% [uu_R2 , omega2_R ] = fminunc(R_q2_p,uu0_R);
% 
% u_R2 = ([1; uu_R2].'*B_Ker_u1.').' ;
% u_R2 = u_R2.*(1/u_R2(1));    % second mode normalize at the first entry
% 
% u_R3 = null([u_R2.' ; u_R1.']);
% u_R3 = u_R3.*(1/u_R3(1));
% 
% U_Rail_p     = [u_R1 , u_R2 , u_R3];
% omega_Rail_p = [R_q_p(u_R1(2:3)) R_q_p(u_R2(2:3)) R_q_p(u_R3(2:3))].^(1/2);
% 
% %pp_raileight
% 
% [U_Rail_p U_eig_p]

%% MODAL ANALYSIS - RAILEIGHT - the new hope (direct method)
syms ur1 ur2

R_q_f  = ( [1; ur1; ur2].' * K_f * [1; ur1; ur2] )  /  ( [1; ur1; ur2].' * M_f * [1; ur1; ur2] ) ;

D1R_q_f = diff(R_q_f,ur1);
D2R_q_f = diff(R_q_f,ur2);

assume(ur1, 'real');
assume(ur2, 'real');
pocheBriciole = vpasolve([ D1R_q_f == 0, D2R_q_f==0],[ur1 ur2]);

U_Rail_f = double([ [1 pocheBriciole.ur1(1) pocheBriciole.ur2(1) ].' ...
                    [1 pocheBriciole.ur1(2) pocheBriciole.ur2(2) ].' ...
                    [1 pocheBriciole.ur1(3) pocheBriciole.ur2(3) ].']);

fun_R_q_f = @(aaa,bbb) ( [1; aaa; bbb].' * K_f * [1; aaa; bbb] )  /  ( [1; aaa; bbb].' * M_f * [1; aaa; bbb] ) ;

for i = 1:3
    wr_f(i)     = sqrt(fun_R_q_f(U_Rail_f(2,i),U_Rail_f(3,i))      )  ;        
end

%% MODAL ANALYSIS - RAILEIGHT - the new hope (direct method)
syms ur1 ur2

R_q_p  = ( [1; ur1; ur2].' * K_p * [1; ur1; ur2] )  /  ( [1; ur1; ur2].' * M_p * [1; ur1; ur2] ) ;

D1R_q_p = diff(R_q_p,ur1);
D2R_q_p = diff(R_q_p,ur2);

assume(ur1, 'real');
assume(ur2, 'real');
pocheBriciole = vpasolve([ D1R_q_p == 0, D2R_q_p==0],[ur1 ur2]);

U_Rail_p = double([ [1 pocheBriciole.ur1(1) pocheBriciole.ur2(1) ].' ...
                    [1 pocheBriciole.ur1(2) pocheBriciole.ur2(2) ].' ...
                    [1 pocheBriciole.ur1(3) pocheBriciole.ur2(3) ].'])
                
fun_R_q_p = @(aaa,bbb) ( [1; aaa; bbb].' * K_p * [1; aaa; bbb] )  /  ( [1; aaa; bbb].' * M_p * [1; aaa; bbb] ) ;

for i = 1:3
    wr_p(i)     = sqrt(fun_R_q_p(U_Rail_p(2,i),U_Rail_p(3,i))      )  ;        
end

pp_raileight

%% MODAL ANALYSIS - Matrix iteration method - full

%Essentially it is a way to find the largest eigen value and the
%corresponding eigenvector. so the procedur eis to find one by one alla the
%omegas. combining MIM-deflation-MIM-deflation-MIM

    %deflation: now we will move this eigenvalue to zero (theoretically), also
    %if it is not zero no problem because this operatin is for remove it from
    %the pole position of the eigenvalues

    %initialization and frst guess
    
    %THE MATRIX DEFLATION DOES NOT AFFECT THE EIGENVECTORS BUT REDUCE THE
    %CHOSEN EIGENVLAUE!
uMIM      = [1 1 1].';
D_0       = K_f\M_f;
U_MIM     = zeros(3);

eig_D_MIM = zeros(3,1);

D_f = D_0;
for i=1:3
%Iterate for the 3 modes
    for j=1:15
        if (j==1)
            uMIM = ones(3,1); %first guess again
        end
        uMIM = D_f*uMIM;
        uMIM = uMIM./uMIM(1);
    end
    %compute i-th eigenvalue (SEE PAGE STANDARD EIGENVALUE PROBLEM ON YOUR NOTEBOOK TO PUT COMPUTATIONS ON REPORT)
    eig_D_MIM(i) = (uMIM.' * D_f * uMIM) / (uMIM.'*uMIM);
    U_MIM(:,i)   = uMIM; 
    %Matrix deflation:
        normalizator = uMIM.' * M_f * uMIM;     % compute the normalization factor
        uMIM = uMIM./( (normalizator)^(1/2) ); % use the normalizaed-to-mass uMIM vector to deflate the matrix D
    D_f = D_f - eig_D_MIM(i)*uMIM*uMIM.'*M_f;
    
    
end

U_MIM_f = U_MIM
[U_MIM_f U_eig_f]

wm_f = sqrt(1./eig_D_MIM)


%%  MODAL ANALYSIS - Matrix iteration method - prop 

%Essentially it is a way to find the largest eigen value and the
%corresponding eigenvector. so the procedur eis to find one by one alla the
%omegas. combining MIM-deflation-MIM-deflation-MIM

    %deflation: now we will move this eigenvalue to zero (theoretically), also
    %if it is not zero no problem because this operatin is for remove it from
    %the pole position of the eigenvalues

    %initialization and frst guess
uMIM      = [1 1 1].';
D_0       = K_p\M_p;
U_MIM     = zeros(3);

eig_D_MIM = zeros(3,1);

D_p = D_0;
for i=1:3
%Iterate for the 3 modes
    for j=1:15
        if (j==1)
            uMIM = ones(3,1); %first guess again
        end
        uMIM = D_p*uMIM;
        uMIM = uMIM./uMIM(1);
    end
    %compute i-th eigenvalue (SEE PAGE STANDARD EIGENVALUE PROBLEM ON YOUR NOTEBOOK TO PUT COMPUTATIONS ON REPORT)
    eig_D_MIM(i) = (uMIM.' * D_p * uMIM) / (uMIM.'*uMIM);
    U_MIM(:,i)   = uMIM; 
    %Matrix deflation:
        normalizator = uMIM.' * M_p * uMIM;     % compute the normalization factor
        uMIM = uMIM./( (normalizator)^(1/2) ); % use the normalizaed-to-mass uMIM vector to deflate the matrix D
    D_p = D_p - eig_D_MIM(i)*uMIM*uMIM.'*M_p;
end

U_MIM_p = U_MIM
[U_MIM_p U_eig_p]

wm_p = sqrt(1./eig_D_MIM)

pp_mim
%% MODAL ANALYSIS - Laplace

%tranfer function with force (voltage-to-force coefficient compensation)
LaPlace_full = (1/g_v_est_f)*tf(sys_guess_full);

sys_guess_full.inputName = 'v';
sys_guess_full.outputName = {'x_1','x_2','x_3'};

freq_range_bode = logspace(0,2,400);



%tranfer function with force (voltage-to-force coefficient compensation)
LaPlace_prop = (1/g_v_est_p)*tf(sys_guess_prop);

pp_laplace

%% MODAL ANALYSIS - analytical proportional damping

% Symbolic matrices

syms m_1 m_2 m_3 k_1 k_2 k_3 c_a c_b;

syms U11 U21 U31 ;
syms U12 U22 U32 ;
syms U13 U23 U33 ;


U_s  = [ U11 U21 U31 ;
         U12 U22 U32 ;
         U13 U23 U33 ];

M_s = [m_1 0 0;  
     0 m_2 0; 
     0 0 m_3];
 
 %K definition fixed
 
 K_s = [k_1   -k_1        0;  
     -k_1 k_1+k_2      -k_2; 
      0    -k_2     k_2+k_3];
  
 C_s = c_a*M_s + c_b*K_s; %called C but it isn't a problem in order of assignments
 
M_s = U_s.'*M_s*U_s;
K_s = U_s.'*K_s*U_s;
C_s = U_s.'*C_s*U_s;



pp_MA_proportionaldamping

%% SINE SWEEP DATA - slow

%fast

Ts_ssf         = tssf(1000) - tssf(999);
data_ssf       = iddata( [ x1ssf  ] , vssf, Ts_ssf )
frosta1 = tfest(data_ssf,[ 6 ],[ 4 ]);

Ts_sss         = tsss(1000) - tsss(999);
data_sss       = iddata( [ x1sss  ] , vsss, Ts_sss )
rispensa1 = tfest(data_sss,[6 ],[ 4 ]);

Ts_ssf         = tssf(1000) - tssf(999);
data_ssf       = iddata( [ x2ssf  ] , vssf, Ts_ssf )
frosta2 = tfest(data_ssf,[ 6 ],[ 3 ]);

Ts_sss         = tsss(1000) - tsss(999);
data_sss       = iddata( [ x2sss  ] , vsss, Ts_sss )
rispensa2 = tfest(data_sss,[6 ],[ 3 ]);

Ts_ssf         = tssf(1000) - tssf(999);
data_ssf       = iddata( [ x3ssf  ] , vssf, Ts_ssf )
frosta3 = tfest(data_ssf,[ 6 ],[ 2 ]);

Ts_sss         = tsss(1000) - tsss(999);
data_sss       = iddata( [ x3sss  ] , vsss, Ts_sss )
rispensa3 = tfest(data_sss,[6 ],[ 2 ]);


%pp_sweep

%% EXTRAS

% frequency response identification with impulse

resp_dataPulse = fft(iddata([ x1, x2, x3 ],v, Ts))

rispensaPulse = tfest(resp_dataPulse,[ 6 6 6 ].', [4 3 2].')

%Friction estimation

close all



