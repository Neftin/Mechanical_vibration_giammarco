
clear all;

%symbolic computations for static gain :

syms k1 k2 k3 gv

K = [k1   -k1        0;
     -k1 k1+k2     -k2;
     0    -k2   k2+k3];
 
g_dc = K\[gv; 0; 0]

%multiply times k_3 supposing is the correct value

g_dc_k3 = expand(simplify(g_dc*k3))

%%  

%loading entire dataset

load('dataMass');

%2.1 estimation of the ratio F/v

%plot the step response 4-times

%see on the plot directly the values, VERY QUICKLY WRITE IT ON REPORT,

%remove t from plot

%ultimate gain

gain_x = (0.0706/16000); %encoder!!!

gain_v = 5.250;

gain_tot = (gain_x/gain_v);

Measured_coeff = zeros(4,3);

%compute the steady state values

aver_sa = 400;

%picking four times averaging the last "average_samples" number of samples

for i = 0:3

Measured_coeff(i+1,1) = ( mean(x1s((950 + 2000*i - aver_sa):(950 + 2000*i))  )) * (gain_x);

Measured_coeff(i+1,2) = ( mean(x2s((950 + 2000*i - aver_sa):(950 + 2000*i))  )) * (gain_x);

Measured_coeff(i+1,3) = ( mean(x3s((950 + 2000*i - aver_sa):(950 + 2000*i))  )) * (gain_x);

end

clear i;

Measured_coeff_means = mean(Measured_coeff);

%first term for the equation k3*x/vstep = g_dc
k3 = 396;
k2 = 770;
k1 = 774;

coeff_for_eq = Measured_coeff_means.*(k3/vstep(950));

g_v = coeff_for_eq(3)

%manually solve the other two equations

%second row of eq
syms r32
R32_exp = double(solve(coeff_for_eq(2) == g_v + g_v*r32 ))

%third row of eq
syms r31
R31_exp = double(solve(coeff_for_eq(1) == g_v + g_v*r31 + g_v*R32_exp ))

R32_nom = k3/k2;
R31_nom = k3/k1;

R31_per = abs((R31_exp - R31_nom)*100 / (R31_exp));
R32_per = abs((R32_exp - R31_exp)*100 / (R31_exp));

g_v_per = abs(100*(g_v - gain_v)/(gain_v));

pp_step

