%post processing estimation
close all;

%% FULL

if ~isempty(figure(3))
     close(3)
end 

 figure(3)
    h = plot(t,v)
    grid on; 
    set(gca,'FontSize',PD/2);
    set(h,'linewidth',1);
    set(h,'color','black')
    xlabel('time [s]','interpreter','latex');
    ylabel('voltage (input) $v$ [V]','interpreter','latex')
    axis([-inf inf -0.2 max(v)*1.2 ]);
    
    axes('position',[0.84 0.35 0.13 0.6])
    h = plot(t,v)
    grid minor; 
    set(h,'linewidth',1);
    set(h,'color','black')
    axis([9.89 10.11 -0.2 max(v)*1.2 ]);
    
    
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 PD PD/3.3]);
    print('report\img\pulse_input','-depsc','-cmyk');%legenda: _per è l'errore percentuale!
%obtained parameters

print2file(m1_f,'report\result\','%3.3f')
print2file(m2_f,'report\result\','%3.3f')
print2file(m3_f,'report\result\','%3.3f')

print2file(c1_f,'report\result\','%3.3f')
print2file(c2_f,'report\result\','%3.3f')
print2file(c3_f,'report\result\','%3.3f')
print2file(c12_f,'report\result\','%1.5f')
print2file(c23_f,'report\result\','%1.5f')

print2file(g_v_est_f,'report\result\','%3.3f')


[ ~, fit_f, ~ ] = compare(comparison_data,sys_guess_full);

y_plot_f      = lsim(sys_guess_full,v,t);
t_plot_f      = (1:length(comparison_data.OutputData) ).*Ts;
yreal_plot_f  = comparison_data.OutputData;
y_resi_plot_f = comparison_data.OutputData - y_plot_f;
for i = 1:3
    figure(4)
    subplot(3,1,i)
    h = plot(t_plot_f,[ yreal_plot_f(:,i) , y_plot_f(:,i), y_resi_plot_f(:,1) ] )
    set(h(1),'color',[0 .3 1]);
    set(h(2),'color',[1 .5 0]);
    set(h(3),'color',[.7 .7 .7]);
%    set(h(3),'LineStyle','-.');
    grid minor;
        xlabel('time [s]','interpreter','latex');
        ylabel(['$x_' num2str(i,1) '$ [m]'],'interpreter','latex');
   hold on;
%     h = area(t_plot_f,y_resi_plot_f(:,i),'LineStyle',':')
%         xlabel('time [s]','interpreter','latex');
%         ylabel(['$x_{' num2str(i,1) 'system}-x_{' num2str(i,1) 'model}$ [m]'],'interpreter','latex');
%     h(1).FaceColor = [0 0 0];
%     h(1).FaceAlpha = 0.3;
   hold off;
   legend('System','Model','Residuals')
   
   
        set(gcf,'PaperUnits','centimeters');
        set(gcf,'PaperPosition',[0 0 PD 1.4*PD]); %A4 has sqrt(2) as ratio!
        print([ 'report\img\compare_f' ],'-depsc','-cmyk');%legenda: _per è l'errore percentuale!
end

for i =(1:3)
    print2file(fit_f(i),'report/results','%3.2f','\n','txt',['fit_f' num2str(i,1)])
end

%% proportional case

%obtained parameters

print2file(m1_p,'report\result\','%3.3f')
print2file(m2_p,'report\result\','%3.3f')
print2file(m3_p,'report\result\','%3.3f')

print2file(ca_p,'report\result\','%3.3f')
print2file(cb_p,'report\result\','%1.3d')

print2file(g_v_est_p,'report\result\','%3.3f')


[ ~, fit_p, ~ ] = compare(comparison_data,sys_guess_prop);

y_plot_p      = lsim(sys_guess_prop,v,t);
t_plot_p      = (1:length(comparison_data.OutputData) ).*Ts;
yreal_plot_p  = comparison_data.OutputData;
y_resi_plot_p = comparison_data.OutputData - y_plot_p;
for i = 1:3
    figure(5)
    subplot(3,1,i)
    h = plot(t_plot_p,[ yreal_plot_p(:,i) , y_plot_p(:,i), y_resi_plot_p(:,1) ] )
    set(h(1),'color',[0 .3 1]);
    set(h(2),'color',[1 .5 0]);
    set(h(3),'color',[.7 .7 .7]);
%    set(h(3),'LineStyle','-.');
    grid minor;
        xlabel('time [s]','interpreter','latex');
        ylabel(['$x_' num2str(i,1) '$ [m]'],'interpreter','latex');
   hold on;
%     h = area(t_plot_p,y_resi_plot_p(:,i),'LineStyle',':')
%         xlabel('time [s]','interpreter','latex');
%         ylabel(['$x_{' num2str(i,1) 'system}-x_{' num2str(i,1) 'model}$ [m]'],'interpreter','latex');
%     h(1).FaceColor = [0 0 0];
%     h(1).FaceAlpha = 0.3;
   hold off;
   legend('System','Model','Residuals')
   
   
        set(gcf,'PaperUnits','centimeters');
        set(gcf,'PaperPosition',[0 0 PD 1.4*PD]); %A4 has sqrt(2) as ratio!
        print([ 'report\img\compare_p' ],'-depsc','-cmyk');%legenda: _per è l'errore percentuale!
end

for i =(1:3)
    print2file(fit_p(i),'report/results','%3.2f','\n','txt',['fit_p' num2str(i,1)])
end

%norm of matrices and relative ca and cb

c_a_m = ca_p*norm(M_p);
c_b_k = cb_p*norm(K_p);

print2file(c_a_m,'report\result\','%3.3f')
print2file(c_b_k,'report\result\','%3.3f')

%%

figure(11)
subplot(1,2,1)
    h = plot(t_plot_p,[yreal_plot_p(:,1) , y_plot_p(:,1),] )
    set(h(1),'color',[0 .3 1]);
    set(h(2),'color',[1 .5 0]);
    axis([ 2 8 -0.002 0.002])
        grid minor;
        xlabel('time [s]','interpreter','latex');
        ylabel(['$x_' num2str(1,1) '$ [m]'],'interpreter','latex');
        
subplot(1,2,2)
    h = plot(t_plot_p,[yreal_plot_p(:,1) , y_plot_p(:,1)] )
    set(h(1),'color',[0 .3 1]);
    set(h(2),'color',[1 .5 0]);
    axis([ 32 38 -0.002 0.002])
        grid minor;
        xlabel('time [s]','interpreter','latex');
        ylabel(['$x_' num2str(1,1) '$ [m]'],'interpreter','latex');
        leg = legend('System','Model');
        set(leg,'FontSize',6)
        
        set(gcf,'PaperUnits','centimeters');
        set(gcf,'PaperPosition',[0 0 PD PD/4]); %A4 has sqrt(2) as ratio!
        print([ 'report\img\dry_fric' ],'-depsc','-cmyk');%legenda: _per è l'errore percentuale!
%%
        
for i = 1:3
    for j = 1:3
        print2file(Cc_f(i,j),'report\result\','%3.3f','\n','txt',['Cdamp' num2str(i,1) num2str(j,1)]);
    end
end
