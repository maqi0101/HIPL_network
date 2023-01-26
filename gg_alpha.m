clc;
clear; 
Sf=24;%%%%%%flower的物种数
Sh=24;%%%%%%%herbivore物种数
Sp=24;%%%%%%%%pollinator物种数
h_fh=10;%%%%%%%herbivore对flower捕食的半饱和常数
h_fp=10;%%%%%%pollinator对flower的半饱和常数
h_pf=10;%%%%%%flower对pollinator的半饱和常数
e_fh=0.6;%%%%%flower到捕食者herbivore的转化系数
% beta=e_ph.*alpha;
mu_h=0.05;%%%%%%herbivore的自然死亡率
mu_p=0.01;%%%%%%pollinator的自然死亡率
%%
P_Nant=0.5;%%%捕食网络嵌套概率
P_Nmut=0.5;%%%互惠网络嵌套概率
P_Qant=0.5;%%%互惠网络模块概率
P_Qmut=0.5;%%%捕食网络模块概率
Cant=0.2;%%%捕食网络连接度
Cmut=0.2;%%%互惠网络连接度
%%
% y0=rand(Sp+Sh+Sm,1)*10;%%%列向量，各物种初始值
t0 =0:0.1:5000;
%%
gg=0:0.05:1;%%%%%%%%干扰项的调节系数
muu=0.5:0.05:1.5;%%%%
rep=500;
cv=zeros(length(gg),length(muu),rep);
for kk=1:rep
    r=0.05*ones(Sf,1);%%%%%%%flower的增长率
    delta_f=0.025*ones(Sf,1);%%%%%flower的密度依赖, Sf*1 matrix
    delta_p=0.025*ones(Sp,1);%%%%%pollinator的密度依赖, Sp*1 matrix
    alpha0=unifrnd(0.2,0.3,Sf,Sh);%%%%%%(i,j)捕食者j对植物i的捕食的潜在系数Sp*Sh matrix
    b_pf=unifrnd(0.2,0.3,Sp,Sf);%%%%%%植物从授粉获益系数
    c_fp=unifrnd(0.2,0.3,Sf,Sp);%%%%%%pollinator从授粉获益系数
    %beta=e_ph.*alpha;
    matrix_FH=get_matrix(Sf,Sh,P_Nant,P_Qant,Cant);%%%捕食关系矩阵
    matrix_FP=get_matrix(Sf,Sp,P_Nmut,P_Qmut,Cmut);%%%互惠关系矩阵
    y0=rand(Sf+Sh+Sp,1)*10;%%%列向量，各物种初始值
    for jj=1:length(muu)
        mu=muu(jj);
        alpha=mu*alpha0;
        beta=e_fh*alpha;
        for ii=1:length(gg)
            g=gg(ii);
            [t,y]=ode45(@(t,y) FHP(t,y,Sf,Sh,Sp,g,h_fh,h_fp,h_pf,r,delta_f,delta_p,alpha,beta,b_pf,c_fp,matrix_FH,matrix_FP,mu_h,mu_p),t0,y0);
            cv(ii,jj,kk)=std(sum(y(end-1000:end,:),2))/mean(sum(y(end-1000:end,:),2));
            fprintf('\nrepeat=%d, row=%d,,list=%d\n\n', kk, ii, jj);
        end
    end
end
%%
stab=zeros(length(gg),length(muu),rep);%%和cv一样规模的零矩阵,多次模拟
for kk=1:rep
    for jj=1:length(muu)
        for ii=1:length(gg)
            if cv(ii,jj,kk)<0.0001  %判断稳定的条件
                stab(ii,jj,kk)=1;
            end
        end
    end
end
P_stab=sum(stab,3)/rep;%稳定的概率
%
figure(1)
imagesc([0.5 1.5],[0 1],P_stab);
axis xy
xlabel('coefficient of antagonistic strength')
ylabel('g')
title('Probability of Stability')
colorbar