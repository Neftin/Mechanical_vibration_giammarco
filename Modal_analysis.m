
%eigenvalues problem!

%Generic Damping case

m1 = Parameters(1);
m2 = Parameters(2);
m3 = Parameters(3);

c1 = Parameters(4);
c2 = Parameters(5);
c3 = Parameters(6);

c12 = Parameters(7);
c23 = Parameters(8);

I = eye(3);

Z = zeros(3);

M = [m1 0 0;  
     0 m2 0; 
     0 0 m3];
 
 %K definition fixed
 
 k1=774;k2=770;k3=396;
 
 K =     [k1   -k1        0;  
         -k1 k1+k2      -k2; 
          0    -k2     k2+k3];
 
  
 C =   [c1+c12       ,-c12                 0;  
       -c12            c2+c23+c12    ,-c23; 
        0            ,-c23             c23+c3];
    
    
syms x;

Omegas = vpasolve( det(K - x*M ) == 0, x);

omegas = double ( (Omegas).^(1/2) )

clear x Omegas;

% rayleight quotient

R_q = @(u) (  transpose( [1 ;u(1); u(2)] ) * K * [1; u(1) ;u(2)] )  /  (  transpose([1; u(1) ;u(2)]) * M * [1; u(1); u(2)]  );

u0_R = [1,2];

[u_R , omega1_R ] = fminunc(R_q,u0_R);

clear u0_R;

u_R = [1 ;u_R(1); u_R(2)]

omega1_R = omega1_R^(1/2)











