
%try to impement the lsqnonlin within the simulation for the output!!!

load('dataMass');

rng default 

rg = (1:length(t));

%data processing and homologation, and choice of range in samples

INPUT_FOR_ESTIMATION = timeseries(v(rg)*gain_v_NEW,t(rg));

OUTPUT_TO_TEND       = timeseries([x1(rg)*gain_x , x2(rg)*gain_x, x3(rg)*gain_x],t(rg));

comparison_data      = iddata([x1(rg)*gain_x , x2(rg)*gain_x , x3(rg)*gain_x], v(rg)*gain_v_NEW,Ts);

%initial guess

global A; %to update by function
global B;
global C;
global D;


%K = [k1 k2 k3];

%P = [ A_mean B1(3:7) B2(4:7) B3(5:7) ]; %first guess by step response estimation (better if manually)

%the zeroes must be dropped!

%first guess

P = [1.5 1.5 1.5 5 5 5 5 5];

%P = [m1,m2,m3,c1,c2,c3,c12,c23]; %MANAGE LOWERBOUNDs AND UPPERBOUNDS, k
%given

lb = [0 0 0 0 0 0 0 0 ];
ub = [2.1 2.1 2.1 10 10 10 10 10];

%%

options = optimoptions('lsqnonlin','Display','iter','StepTolerance',0.000000001,'OptimalityTolerance',0.000001);

Parameters = lsqnonlin(@errors_direct,P,lb,ub, options) %the core of estimation

%%

sys_guess = ss(A,B,C,D);

figure(1);
compare(comparison_data,sys_guess)


