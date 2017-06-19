%pp_laplace

freq_range_bode = logspace(0,2,400);

close all;

figure(5);
bode((1/(g_v_est_f))*LaPlace_full,freq_range_bode);
hold on;
grid minor;

figure(5);
bode((1/(g_v_est_f))*LaPlace_prop,freq_range_bode);
leg = legend('Free damping','proportional damping')
set(leg,'FontSize',6)
hold off;
grid minor;



    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 PD 1.4*PD]);
    print('report\img\bode_fp','-depsc','-cmyk');
    
