%loading entire dataset

load('dataMass');

% 2.1 estimation of the ratio F/v

% plot the step response 4-times

%see on the plot directly the values, VERY QUICKLY WRITE IT ON REPORT,

%remove t from plot

%ultimate gain

gain_x = (0.0706/16000); %encoders??!

gain_v = 5.250;

gain_tot = (gain_x/gain_v);

Measured_coeff = zeros(4,3);

for i = 0:3

Measured_coeff(i+1,1) = ( x1s(950 + 2000*i) /vstep(950) ) * (gain_tot);

Measured_coeff(i+1,2) = ( x2s(950 + 2000*i) /vstep(950) ) * (gain_tot);

Measured_coeff(i+1,3) = ( x3s(950 + 2000*i) /vstep(950) ) * (gain_tot);

end

clear i;

Measured_coeff_means = mean(Measured_coeff);

