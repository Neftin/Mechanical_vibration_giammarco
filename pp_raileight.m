close all;
%pp_Raileight

%free damping

%Raileight quotient to plot
R_q_peova_f = @(x,y) ( [1; x; y].' * K_f * [1; x; y] )  /  ( [1; x; y].' * M_f * [1; x; y] ) ;

% R_q_peova_f = @(xx) ( [1; xx(1); xx(2)].' * K_f * [1; xx(1); xx(2)] )  /  ( [1; xx(1); xx(2)].' * M_f * [1; xx(1); xx(2)] ) ;
figure(6);

fcontour(R_q_peova_f,[-3 3 -2.5 2.5],'LevelList',[logspace(-3,6,300)]); hold on;
grid minor;
    xlabel('$\alpha$','interpreter','latex');
    ylabel('$\beta$','interpreter','latex')
figure(6); h = plot(U_Rail_f(2,1),U_Rail_f(3,1),'+',U_Rail_f(2,2),U_Rail_f(3,2),'+',U_Rail_f(2,3),U_Rail_f(3,3),'+')
set(h,'MarkerSize',14);
set(h,'linewidth',2);
hold off;

    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 PD PD/1.7]);
    print('report\img\contour_f','-depsc','-cmyk');

%proportional

%Raileight quotient to plot
R_q_peova_p = @(x,y) ( [1; x; y].' * K_p * [1; x; y] )  /  ( [1; x; y].' * M_p * [1; x; y] ) ;

% R_q_peova_p = @(xx) ( [1; xx(1); xx(2)].' * K_p * [1; xx(1); xx(2)] )  /  ( [1; xx(1); xx(2)].' * M_p * [1; xx(1); xx(2)] ) ;
figure(7);

fcontour(R_q_peova_p,[-3 3 -2.5 2.5],'LevelList',[logspace(-3,6,300)]); hold on;
grid minor;
    xlabel('$\alpha$','interpreter','latex');
    ylabel('$\beta$','interpreter','latex')
figure(7); h = plot(U_Rail_p(2,1),U_Rail_p(3,1),'+',U_Rail_p(2,2),U_Rail_p(3,2),'+',U_Rail_p(2,3),U_Rail_p(3,3),'+')
set(h,'MarkerSize',14);
set(h,'linewidth',2);

hold off;

    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 PD PD/1.7]);
    print('report\img\contour_p','-depsc','-cmyk');
    
    
for i = 1:3
    print2file(wr_f(i),'report\result\','%3.4f','\n','txt',['wr_f_' num2str(i,1)])
        for j = 1:3
            print2file(U_Rail_f(j,i),'report\result\','%3.4f','\n','txt',['Ur_f_' num2str(j,1) num2str(i,1)])
        end
end




for i = 1:3
    print2file(wr_p(i),'report\result\','%3.4f','\n','txt',['wr_p_' num2str(i,1)])
            for j = 1:3
            print2file(U_Rail_p(j,i),'report\result\','%3.4f','\n','txt',['Ur_p_' num2str(j,1) num2str(i,1)])
            end
end