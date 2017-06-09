



% rayleight quotient

[1 u(1) u(2)]

R_q = @u(  transpose([1 u(1) u(2)]) * K * [1 u(1) u(2)] )  /  ( transpose([1 u(1) u(2)]) * M * [1 u(1) u(2)]  );

u0_R = [1,2] %change it, is just an example

[u_R , omega1_R ] = fminunc(R_q,u0_R)
