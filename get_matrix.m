function matrix_FA = get_matrix(SSf,SSa,P_N,P_Q,C0)
f_pl=(1:SSf).^-2;%%power-law distribution of degree-2
a_pl=(1:SSa).^-2;
%
A_site=randsrc(1,SSa,[1:SSf;f_pl/sum(f_pl)]);%%%��������������ÿ�������з���һ��ֵ�������������������������໥���õĸ���
F_site=randsrc(1,SSf,[1:SSa;a_pl/sum(a_pl)]);%%%��������������ÿ��ֲ���з���һ��ֵ�������������������������໥���õĸ���
%%
matrix_FA=zeros(SSf,SSa);
%
Ln_total=0;%%%�����������ʼ��������
while Ln_total<round(C0*SSf*SSa)
    Type=randi(2);   %%% ѡ���һ���������ڵ��л�,1:P->A; 2:A->P
    if Type==1%��ѡֲ����������
        if rand<=P_N%%%ֲ���з���Ƕ�׵ĸ���
            f_si=randsrc(1,1,[1:SSf;F_site/sum(F_site)]);%%ѡȡһ��ֲ��
        else
            f_si=randsrc(1,1,[1:SSf;1/SSf*ones(1,SSf)]);
        end
        si_part=ceil(4*f_si/SSf);%%%ȷ����ֲ�����ڵ�ģ�飨1,2,3,4)
        if rand<=P_Q%%%%����ģ�黯�ĸ���
            A_rang=1+(si_part-1)*SSa/4:si_part*SSa/4;%%%ȷ���۷�ķ�Χ
            A_dens=A_site(1+(si_part-1)*SSa/4:si_part*SSa/4);%%%�÷�Χ�µ���ֵ
        else
            A_rang=1:SSa;
            A_dens=A_site;
        end
        if rand<=P_N%%%%���﷢��Ƕ�׵ĸ���
            a_si=randsrc(1,1,[A_rang;A_dens/sum(A_dens)]);%%%ȷ��һ���۷�
        else
            a_si=randsrc(1,1,[A_rang;1/length(A_rang)*ones(1,length(A_rang))]);
        end
        matrix_FA(f_si,a_si)=1;%%%��ֲ��Ͷ�������
    else%��ѡ������ѡֲ��
        if rand<=P_N%%%�����з���Ƕ�׵ĸ���
            a_si=randsrc(1,1,[1:SSa;A_site/sum(A_site)]);%%ѡȡһ������
        else
            a_si=randsrc(1,1,[1:SSa;1/SSa*ones(1,SSa)]);
        end
        si_part=ceil(4*a_si/SSa);%%%ȷ��ֲ�����ڵ�ģ��
        if rand<=P_Q%%%%����ģ��ĸ���
            P_rang=1+(si_part-1)*SSf/4:si_part*SSf/4;%%%ȷ��ֲ��ķ�Χ
            P_dens=F_site(1+(si_part-1)*SSf/4:si_part*SSf/4);%%%�÷�Χ�µ���ֵ
        else
            P_rang=1:SSf;
            P_dens=F_site;
        end
        if rand<=P_N%%%%ֲ�﷢��Ƕ�׵ĸ���
            f_si=randsrc(1,1,[P_rang;P_dens/sum(P_dens)]);%%%ȷ��һ��ֲ��
        else
            f_si=randsrc(1,1,[P_rang;1/length(P_rang)*ones(1,length(P_rang))]);
        end
        matrix_FA(f_si,a_si)=1;%%%�������ֲ������
    end
    Ln_total=sum(sum(matrix_FA));%%%���������ڸ�ʱ�̵���������
end
