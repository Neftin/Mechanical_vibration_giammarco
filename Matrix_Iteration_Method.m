%% Matrix Iteration Method

% Paramentri stimati
m1 = Parameters(1);
m2 = Parameters(2);
m3 = Parameters(3);

c1 = Parameters(4);
c2 = Parameters(5);
c3 = Parameters(6);

c12 = Parameters(7);
c23 = Parameters(8);
%
k1 = 770;
k2 = 774;
k3 = 396;
%
M = [m1 0  0;  
     0  m2 0; 
     0  0  m3];
%
K = [ k1    -k1       0;  
      -k1   k1+k2     -k2; 
       0     -k2     k2+k3];
 %
 % Iterations
 n = 10 ;
 %
 %% 1st Mode
 % FIRST Deflated Matrix
 D1=K\M;
 % Empty arrays and matrices to store results in.
 X1         = zeros(3,n);
 lambda1    = zeros(1,n);
 omega1     = zeros(1,n);
 % First guess
 X1(:,1)    = [1; 1; 1];
 lambda1(1) = 1;
 omega1(1)  = 1; 
 %
 % Cycle that evaluate  the FIRST mode and the FIRST natural frequency
 for i = 1 : n-1
    X1(:,i+1)  = D1*X1(:,i); 
    lambda1(i+1) = X1(1,i+1);
    omega1(i+1)  = sqrt(1/lambda1(i+1));
    X1(:,i+1)  = X1(:,i+1)/lambda1(i+1);
 end
 error = omega1(n) - omega1(n-1);
 
 %% 2nd Mode
 syms alpha1;
 tmp1 = alpha1.*X1(:,n);
 Alpha1 = vpasolve( transpose(tmp1)*M*tmp1 == 1, alpha1);
 Alpha1 = double(Alpha1(2,1));
 %
 % SECOND Deflated matrix
 D2 = D1 - lambda1(n)*Alpha1*Alpha1*X1(:,n)*transpose(X1(:,n))*M;
  % Empty arrays and matrices to store results in.
 X2         = zeros(3,n);
 lambda2    = zeros(1,n);
 omega2     = zeros(1,n);
 % First guess
 X2(:,1)    = [1; 1; 1];
 lambda2(1) = 1;
 omega2(1)  = 1; 
  %
 % Cycle that evaluate  the FIRST mode and the FIRST natural frequency
 for i = 1 : n-1
    X2(:,i+1)  = D2*X2(:,i); 
    lambda2(i+1) = X2(1,i+1);
    omega2(i+1)  = sqrt(1/lambda2(i+1));
    X2(:,i+1)  = X2(:,i+1)/lambda2(i+1);
 end
  %% 3nd Mode
 syms alpha2;
 tmp2 = alpha2.*X2(:,n);
 Alpha2 = vpasolve( transpose(tmp2)*M*tmp2 == 1, alpha2);
 Alpha2 = double(Alpha2(2,1));
 %
 % SECOND Deflated matrix
 D3 = D2 - lambda2(n)*Alpha2*Alpha2*X2(:,n)*transpose(X2(:,n))*M;
  % Empty arrays and matrices to store results in.
 X3         = zeros(3,n);
 lambda3    = zeros(1,n);
 omega3     = zeros(1,n);
 % First guess
 X3(:,1)    = [1; 1; 1];
 lambda3(1) = 1;
 omega3(1)  = 1; 
  %
 % Cycle that evaluate  the FIRST mode and the FIRST natural frequency
 for i = 1 : n-1
    X3(:,i+1)  = D3*X3(:,i); 
    lambda3(i+1) = X3(1,i+1);
    omega3(i+1)  = sqrt(1/lambda3(i+1));
    X3(:,i+1)  = X3(:,i+1)/lambda3(i+1);
 end
 
 %% Resluts
 U1 = X1(:,n);
 U2 = X2(:,n);
 U3 = X3(:,n);
 U  = [U1, U2, U3];
 %
 formatted_string1 = sprintf(' The RESONANCE FREQUENCIES are:\n\nomega1 =  %0.5f\nomega2 =  %0.5f\nomega3 =  %0.5f\n', omega1(n),omega2(n),omega3(n));
 disp(formatted_string1);
 formatted_string2 = sprintf(' The  MODE Matrix is:\n');
 disp(formatted_string2);
 disp(U);
 