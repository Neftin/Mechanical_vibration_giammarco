%post processing estimation


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
    h = plot(t_plot_f,[ yreal_plot_f(:,i) , y_plot_f(:,i)  ] )
    set(h(1),'color',[0 .3 1]);
    set(h(2),'color',[1 .5 0]);
    grid minor;
        xlabel('time [s]','interpreter','latex');
        ylabel(['$x_' num2str(i,1) '$ [m]'],'interpreter','latex');
   hold on;
    h = area(t_plot_f,y_resi_plot_f(:,i),'LineStyle',':')
        xlabel('time [s]','interpreter','latex');
        ylabel(['$x_{' num2str(i,1) 'system}-x_{' num2str(i,1) 'model}$ [m]'],'interpreter','latex');
    h(1).FaceColor = [0 0 0];
    h(1).FaceAlpha = 0.3;
   hold off;
        set(gcf,'PaperUnits','centimeters');
        set(gcf,'PaperPosition',[0 0 PD 1.4*PD]); %A4 has sqrt(2) as ratio!
        print([ 'report\img\compare_' num2str(i,1) ],'-depsc','-cmyk');%legenda: _per è l'errore percentuale!
end
%plot (comparisons)

%%

