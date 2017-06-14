%post processing for step response

PD    = 18;
RATIO = 2;% for a fortunate coincidence, a page contains ratios figures

M_s_matrix_tex = latex(M_s);

K_s_matrix_tex = {'1','2','3'};

for i = 1:3
    K_s_matrix_tex{i} = latex(K_s(:,i));
    C_s_matrix_tex{i} = latex(K_s(:,i));
end

for i = ['1','2','3']
    for j = ['1','2','3']
        M_s_matrix_tex= strrep(M_s_matrix_tex,['u' i j ],['u_{' i j '}']);
        K_s_matrix_tex= strrep(K_s_matrix_tex,['u' i j ],['u_{' i j '}']);
        C_s_matrix_tex= strrep(C_s_matrix_tex,['u' i j ],['u_{' i j '}']);
    end
end


print2file(M_s_matrix_tex,'report\result\','%s','\n')

for i = 1:3
   print2file(K_s_matrix_tex{i},'report\result\','%s','\n','txt',['K_s_matrix_tex' num2str(i,1)]);
   print2file(C_s_matrix_tex{i},'report\result\','%s','\n','txt',['C_s_matrix_tex' num2str(i,1)]);
end

