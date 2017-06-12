
%pp_Raileight

R_q_peova = @(x,y) ( [1; x; y].' * K * [1; x; y] )  /  ( [1; x; y].' * M * [1; x; y] ) ;

figure(11); fcontour(R_q_peova,[-2.5 2.5 -2.5 2.5],'LevelList',[logspace(-3,6,300)]); hold on;
%mancano assi nomi e significati
figure(11); h = plot(u_R1(2),u_R1(3),'+',u_R2(2),u_R2(3),'+',u_R3(2),u_R3(3),'+')
set(h,'MarkerSize',20);
set(h,'linewidth',2);
hold off;

figure(12)
fplot(R_q2,[-1,1])
