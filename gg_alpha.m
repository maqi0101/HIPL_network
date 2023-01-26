clc;
clear; 
Sf=24;%%%%%%flower��������
Sh=24;%%%%%%%herbivore������
Sp=24;%%%%%%%%pollinator������
h_fh=10;%%%%%%%herbivore��flower��ʳ�İ뱥�ͳ���
h_fp=10;%%%%%%pollinator��flower�İ뱥�ͳ���
h_pf=10;%%%%%%flower��pollinator�İ뱥�ͳ���
e_fh=0.6;%%%%%flower����ʳ��herbivore��ת��ϵ��
% beta=e_ph.*alpha;
mu_h=0.05;%%%%%%herbivore����Ȼ������
mu_p=0.01;%%%%%%pollinator����Ȼ������
%%
P_Nant=0.5;%%%��ʳ����Ƕ�׸���
P_Nmut=0.5;%%%��������Ƕ�׸���
P_Qant=0.5;%%%��������ģ�����
P_Qmut=0.5;%%%��ʳ����ģ�����
Cant=0.2;%%%��ʳ�������Ӷ�
Cmut=0.2;%%%�����������Ӷ�
%%
% y0=rand(Sp+Sh+Sm,1)*10;%%%�������������ֳ�ʼֵ
t0 =0:0.1:5000;
%%
gg=0:0.05:1;%%%%%%%%������ĵ���ϵ��
muu=0.5:0.05:1.5;%%%%
rep=500;
cv=zeros(length(gg),length(muu),rep);
for kk=1:rep
    r=0.05*ones(Sf,1);%%%%%%%flower��������
    delta_f=0.025*ones(Sf,1);%%%%%flower���ܶ�����, Sf*1 matrix
    delta_p=0.025*ones(Sp,1);%%%%%pollinator���ܶ�����, Sp*1 matrix
    alpha0=unifrnd(0.2,0.3,Sf,Sh);%%%%%%(i,j)��ʳ��j��ֲ��i�Ĳ�ʳ��Ǳ��ϵ��Sp*Sh matrix
    b_pf=unifrnd(0.2,0.3,Sp,Sf);%%%%%%ֲ����ڷۻ���ϵ��
    c_fp=unifrnd(0.2,0.3,Sf,Sp);%%%%%%pollinator���ڷۻ���ϵ��
    %beta=e_ph.*alpha;
    matrix_FH=get_matrix(Sf,Sh,P_Nant,P_Qant,Cant);%%%��ʳ��ϵ����
    matrix_FP=get_matrix(Sf,Sp,P_Nmut,P_Qmut,Cmut);%%%���ݹ�ϵ����
    y0=rand(Sf+Sh+Sp,1)*10;%%%�������������ֳ�ʼֵ
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
stab=zeros(length(gg),length(muu),rep);%%��cvһ����ģ�������,���ģ��
for kk=1:rep
    for jj=1:length(muu)
        for ii=1:length(gg)
            if cv(ii,jj,kk)<0.0001  %�ж��ȶ�������
                stab(ii,jj,kk)=1;
            end
        end
    end
end
P_stab=sum(stab,3)/rep;%�ȶ��ĸ���
%
figure(1)
imagesc([0.5 1.5],[0 1],P_stab);
axis xy
xlabel('coefficient of antagonistic strength')
ylabel('g')
title('Probability of Stability')
colorbar