close all;
%pp_mim

for i = 1:3
    print2file(wm_f(i),'report\result\','%3.4f','\n','txt',['wm_f_' num2str(i,1)])
        for j = 1:3
            print2file(U_MIM_f(j,i),'report\result\','%3.4f','\n','txt',['Um_f_' num2str(j,1) num2str(i,1)])
        end
end




for i = 1:3
    print2file(wm_p(i),'report\result\','%3.4f','\n','txt',['wm_p_' num2str(i,1)])
            for j = 1:3
            print2file(U_MIM_p(j,i),'report\result\','%3.4f','\n','txt',['Um_p_' num2str(j,1) num2str(i,1)])
            end
end