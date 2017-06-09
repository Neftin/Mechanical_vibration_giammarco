function [ F ] = errors_statespace_proportional( P )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

global v gain_v_NEW x1 x2 x3 gain_x t;

m1 = P(1);
m2 = P(2);
m3 = P(3);

ca = P(4);
cb = P(5);

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
     
 A = [Z I; -M\K -M\C]; %left divide for the inverse
 B = [Z; M\I];
 C = [I Z];
 D = Z;
 
ff2 = zeros(7401,1, 'double');
ff3 = zeros(7401,1, 'double');
F = [v*gain_v_NEW ff2 ff3];    %tocca fargliele entrare nella funzione di merda ora vedremo come

W = ones(length(x1),1)*(1/10); %reducing coefficient for error in static friction
active_range_samples = 650; %real active part

for i = 0:3
    
    W((1 + i*1999):((i*1999) +active_range_samples )) = ones(1,active_range_samples);
    
end

X0 = [0 0 0 0 0 0];
 
sys = ss(A, B, C, D);
y = lsim(sys, F, t, X0);

%aggiungi una ponderazione dell'errore tutta sulla parte non piatta!

F = [ ( y(:,1) - x1*gain_x ).*W (y(:,2) - x2*gain_x).*W (y(:,3) - x3*gain_x).*W ];

end

