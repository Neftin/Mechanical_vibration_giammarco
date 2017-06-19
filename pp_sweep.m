
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