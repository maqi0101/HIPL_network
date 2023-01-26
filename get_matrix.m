function matrix_FA = get_matrix(SSf,SSa,P_N,P_Q,C0)
f_pl=(1:SSf).^-2;%%power-law distribution of degree-2
a_pl=(1:SSa).^-2;
%
A_site=randsrc(1,SSa,[1:SSf;f_pl/sum(f_pl)]);%%%给二分子网络中每个动物中分配一个值，用来决定了其在子网络中相互作用的概率
F_site=randsrc(1,SSf,[1:SSa;a_pl/sum(a_pl)]);%%%给二分子网络中每个植物中分配一个值，用来决定了其在子网络中相互作用的概率
%%
matrix_FA=zeros(SSf,SSa);
%
Ln_total=0;%%%二分子网络初始链接总数
while Ln_total<round(C0*SSf*SSa)
    Type=randi(2);   %%% 选择第一个物种所在的行会,1:P->A; 2:A->P
    if Type==1%先选植物再连动物
        if rand<=P_N%%%植物中发生嵌套的概率
            f_si=randsrc(1,1,[1:SSf;F_site/sum(F_site)]);%%选取一个植物
        else
            f_si=randsrc(1,1,[1:SSf;1/SSf*ones(1,SSf)]);
        end
        si_part=ceil(4*f_si/SSf);%%%确定该植物所在的模块（1,2,3,4)
        if rand<=P_Q%%%%发生模块化的概率
            A_rang=1+(si_part-1)*SSa/4:si_part*SSa/4;%%%确定蜜蜂的范围
            A_dens=A_site(1+(si_part-1)*SSa/4:si_part*SSa/4);%%%该范围下的数值
        else
            A_rang=1:SSa;
            A_dens=A_site;
        end
        if rand<=P_N%%%%动物发生嵌套的概率
            a_si=randsrc(1,1,[A_rang;A_dens/sum(A_dens)]);%%%确定一个蜜蜂
        else
            a_si=randsrc(1,1,[A_rang;1/length(A_rang)*ones(1,length(A_rang))]);
        end
        matrix_FA(f_si,a_si)=1;%%%将植物和动物连接
    else%先选动物再选植物
        if rand<=P_N%%%动物中发生嵌套的概率
            a_si=randsrc(1,1,[1:SSa;A_site/sum(A_site)]);%%选取一个动物
        else
            a_si=randsrc(1,1,[1:SSa;1/SSa*ones(1,SSa)]);
        end
        si_part=ceil(4*a_si/SSa);%%%确定植物所在的模块
        if rand<=P_Q%%%%发生模块的概率
            P_rang=1+(si_part-1)*SSf/4:si_part*SSf/4;%%%确定植物的范围
            P_dens=F_site(1+(si_part-1)*SSf/4:si_part*SSf/4);%%%该范围下的数值
        else
            P_rang=1:SSf;
            P_dens=F_site;
        end
        if rand<=P_N%%%%植物发生嵌套的概率
            f_si=randsrc(1,1,[P_rang;P_dens/sum(P_dens)]);%%%确定一个植物
        else
            f_si=randsrc(1,1,[P_rang;1/length(P_rang)*ones(1,length(P_rang))]);
        end
        matrix_FA(f_si,a_si)=1;%%%将动物和植物连接
    end
    Ln_total=sum(sum(matrix_FA));%%%二分网络在该时刻的链接总数
end
