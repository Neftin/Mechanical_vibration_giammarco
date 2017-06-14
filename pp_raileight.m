
%pp_Raileight

%Raileight quotient to plot
R_q_peova_f = @(x,y) ( [1; x; y].' * K_f * [1; x; y] )  /  ( [1; x; y].' * M_f * [1; x; y] ) ;

% R_q_peova_f = @(xx) ( [1; xx(1); xx(2)].' * K_f * [1; xx(1); xx(2)] )  /  ( [1; xx(1); xx(2)].' * M_f * [1; xx(1); xx(2)] ) ;
figure(11); fcontour(R_q_peova_f,[-2.5 2.5 -2.5 2.5],'LevelList',[logspace(-3,6,300)]); hold on;
%mancano assi nomi e significati
% xy0 = [1 1];
% 
% [x1 aaaa] = fminunc(R_q_peova_f,xy0);
% 
% xx_R1 = [1 x1(1) x1(2)];
% 
% H = null(xx_R1);
% 
% Indispensabile = @(yy) (([yy 1] * H.') * K_f * ([yy 1]* H.').') /...
%                        (([yy 1] * H.') * M_f * ([yy 1]* H.').')
% 
%            
% [y1 aaaa] = fminunc(Indispensabile,1);
% 
%  fallimento = [ y1 1]*H.'      
%  
%  fallimento = fallimento./fallimento(1)
        




figure(11); h = plot(u_R1(2),u_R1(3),'+',u_R2(2),u_R2(3),'+',u_R3(2),u_R3(3),'+')
set(h,'MarkerSize',20);
set(h,'linewidth',2);
hold off;

% figure(12);
% fplot(R_q2_f,[-1,1])
