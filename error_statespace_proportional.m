function [ F ] = error_statespace_proportional( P )
%lo spazio di stato è l'ABC dei controlli... usare solo nello script

global v t x1 x2 x3 gain_x;

m1      = P(1);
m2      = P(2);
m3      = P(3);

ca      = P(4);
cb      = P(5);

g_v_est = P(6);





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
    
 
 
 b = [g_v_est 0 0].';
 
 A = [Z I; -M\K -M\C];    %left divide for the inverse
 B = [Z(:,1); M\b]; %single input
 C = [I Z];               %multiple output
 D = Z(:,1);

W = ones(length(x1),1)*(1/10); %reducing coefficient for error in static friction
active_range_samples = 650; %real active part

for i = 0:3
    
    W((1 + i*1999):((i*1999) +active_range_samples )) = ones(1,active_range_samples);
    
end

X0 = [0 0 0 0 0 0];

sys = ss(A, B, C, D);

y = lsim(sys, v, t, X0);

%aggiungi una ponderazione dell'errore tutta sulla parte non piatta!

F = [ ( y(:,1) - x1*gain_x ).*W; (y(:,2) - x2*gain_x).*W; (y(:,3) - x3*gain_x).*W ];

end

