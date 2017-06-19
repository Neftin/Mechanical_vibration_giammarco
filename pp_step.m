close all;
%post processing for step response

PD    = 18;
RATIO = 2;% for a fortunate coincidence, a page contains ratios figures

%11111111111111111
colPl = smartColorPlot(4,200,0.7,'perceived')


for i = 0:3

temp1(i+1) = ( x1s(950 + 2000*i)) * (gain_x);
ttemp1(i+1) = t(950 + 2000*i);
temp2(i+1) = ( x2s(950 + 2000*i)) * (gain_x);

temp3(i+1) = ( x3s(950 + 2000*i)) * (gain_x);

end

figure(1)
h = plot(t,x1s*(gain_x),ttemp1,temp1,'x',...
         t,x2s*(gain_x),ttemp1,temp2,'x',...
         t,x3s*(gain_x),ttemp1,temp3,'x')
    for i = 1:4
        hh = line([ttemp1(i),ttemp1(i)],ylim,'Linestyle',':','color','red','Linewidth',1)
    end
    grid on; 
    set(gca,'FontSize',PD/2);
    set(h,'linewidth',1);
    set(h(2),'MarkerSize',0.5*PD);
    set(h(4),'MarkerSize',0.5*PD);
    set(h(6),'MarkerSize',0.5*PD);
    xlabel('time [s]','interpreter','latex');
    ylabel('positions $x_i$ [m]','interpreter','latex')
    h = legend('x_1(t)','x_1 steady state','x_2(t)','x_2 steady state','x_3(t)','x_3 steady state')
    set(h,'fontsize',5,'fontweight','b') ;
    set(h,'position',[.8,.81,.1,.14])
    set(h,'box','on')
 figure(1)   

     %   printing parameters (last!)
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 PD PD/2]);
    print('report\img\stst1','-depsc','-cmyk');

 figure(2)
    h = plot(t,vstep)
    grid on; 
    set(gca,'FontSize',PD/2);
    set(h,'linewidth',1);
    set(h,'color','black')
    xlabel('time [s]','interpreter','latex');
    ylabel('voltage (input) $v$ [V]','interpreter','latex')
    axis([-inf inf -0.3 0.6 ]);

    
         %   printing parameters (last!)
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 PD PD/3.3]);
    print('report\img\step_input','-depsc','-cmyk');%legenda: _per è l'errore percentuale!
 
results = [R31_exp R31_nom R32_exp R32_nom R31_per R32_per]
names    = {'R31_exp' 'R31_nom' 'R32_exp' 'R32_nom' 'R31_per' 'R32_per'};

%create a file for each result, only to be very fancy
for i = 1:6
    fileID = fopen( ['report\result\' names{i} '.txt'],'wt');
        if (i > 4)
            fprintf(fileID,'%3.1f\n',results(i));
        else
            fprintf(fileID,'%3.4f\n',results(i));
        end
    fclose(fileID);
end

print2file(g_v,'report\result\','%3.4f')
print2file(gain_v,'report\result\','%3.4f')
print2file(g_v_per,'report\result\','%2.1f')


