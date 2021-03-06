%post processing for step response

PD    = 18;
RATIO = 2;% for a fortunate coincidence, a page contains ratios figures

%The nondiagonal azzerator
    for i = 1:3
        for j = 1:3
            if(i ~= j)
                M_s(i,j) = 0;
                K_s(i,j) = 0;
                C_s(i,j) = 0;
            end
        end
    end



M_s_matrix_tex = latex(M_s);

K_s_matrix_tex = {'1','2','3'};

for i = 1:3
    K_s_matrix_tex{i} = latex(K_s(:,i));
    C_s_matrix_tex{i} = latex(C_s(:,i));
end

for i = ['1','2','3']
    for j = ['1','2','3']
        M_s_matrix_tex= strrep(M_s_matrix_tex,['U' i j ],['U_{' i j '}']);
        K_s_matrix_tex= strrep(K_s_matrix_tex,['U' i j ],['U_{' i j '}']);
        C_s_matrix_tex= strrep(C_s_matrix_tex,['U' i j ],['U_{' i j '}']);
    end
end


print2file(M_s_matrix_tex,'report\result\','%s','\n')

for i = 1:3
   print2file(K_s_matrix_tex{i},'report\result\','%s','\n','txt',['K_s_matrix_tex' num2str(i,1)]);
   print2file(C_s_matrix_tex{i},'report\result\','%s','\n','txt',['C_s_matrix_tex' num2str(i,1)]);
end
%% 
zeta_n = (ca_p + cb_p.*w_p.^2)./(2*w_p); %dampingness
omega_n = w_p.*sqrt(1-zeta_n.^2); %

for i = 1:3
    print2file(omega_n(i),'report\result\','%3.3f','\n','txt',['omega_n' num2str(i,1)]);
    print2file(zeta_n(i),'report\result\','%3.3f','\n','txt',['zeta_n' num2str(i,1)]);
end





