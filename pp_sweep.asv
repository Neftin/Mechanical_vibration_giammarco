close all;
%pp_ssweep


%plot sine sweep data
figure(8)
subplot(1,2,1)
zsss = fftshift(fft(vsss));
N = length(zsss);
f = (1/(Ts_sss))*[-N/2:N/2-1]/N;
plot(f,abs(zsss).*Ts_sss)
axis([0 20 0 0.5]);
    xlabel('frequency [Hz]','interpreter','latex')
    ylabel('Amplitude [v]','interpreter','latex')
grid minor;

figure(8)

subplot(1,2,2)
zssf = fftshift(fft(vssf));
N = length(zssf);
f = (1/(Ts_ssf))*[-N/2:N/2-1]/N;
plot(f,abs(zssf).*Ts_ssf)
axis([0 20 0 0.5]);
grid minor;
    xlabel('frequency [Hz]','interpreter','latex')
    ylabel('Amplitude [v]','interpreter','latex')
    
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 PD PD/2.5]);
    print('report\img\sine_sweep_spectra','-depsc','-cmyk');
    
    
 %plot sine sweep data
figure(9)
subplot(1,2,1)
plot(tsss,[x1sss,x2sss,x3sss].*(gain_x))
    xlabel('time [s]','interpreter','latex')
    ylabel('Amplitude [m]','interpreter','latex')
grid minor;
legend('x_1','x_2','x_3')

figure(9)
subplot(1,2,2)
plot(tssf,[x1ssf,x2ssf,x3ssf].*(gain_x))
grid minor;
    xlabel('time [s]','interpreter','latex')
    ylabel('Amplitude [m]','interpreter','latex')
    
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 PD PD/2.5]);
    print('report\img\sine_sweep_time','-depsc','-cmyk');   
 
figure(10);
bode([frosta1 frosta2 frosta3],freq_range_bode);
hold on;
figure(10);
bode([rispensa1 rispensa2 rispensa3],freq_range_bode);
hold off;
grid minor;

    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 1.7*PD 0.7*PD/1.8]);
    print('report\img\ssbode1','-depsc','-cmyk');   
    
figure(11);
bode(frosta2,freq_range_bode);
hold on;
figure(11);
bode(rispensa2,freq_range_bode);
hold off;
grid minor;
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 PD PD/2]);
    print('report\img\ssbode2','-depsc','-cmyk');   
    
figure(12);
bode(frosta3,freq_range_bode);
hold on;
figure(12);
bode(rispensa3,freq_range_bode);
hold off;
grid minor;
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 PD PD/2]);
    print('report\img\ssbode3','-depsc','-cmyk');   