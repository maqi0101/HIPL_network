clc;
clear; 
Sf=24;%%%%%%flower��������
Sh=24;%%%%%%%herbivore������
Sp=24;%%%%%%%%pollinator������
h_fh=10;%%%%%%%herbivore��flower��ʳ�İ뱥�ͳ���
h_fp=10;%%%%%%pollinator��flower�İ뱥�ͳ���
h_pf=10;%%%%%%flower��pollinator�İ뱥�ͳ���
r=0.05*ones(Sf,1);%%%%%%%flower��������
delta_f=0.025*ones(Sf,1);%%%%%flower���ܶ�����, Sf*1 matrix
delta_p=0.025*ones(Sp,1);%%%%%pollinator���ܶ�����, Sf*1 matrix
alpha=unifrnd(0.2,0.3,Sf,Sh);%%%%%%(i,j)herbivore_j��flower_i�Ĳ�ʳϵ��Sf*Sh matrix
b_pf=unifrnd(0.2,0.3,Sp,Sf);%%%%%%flower���ڷۻ���ϵ��
c_fp=unifrnd(0.2,0.3,Sf,Sp);%%%%%%pollinator���ڷۻ���ϵ��
e_fh=0.6;%%%%%flower��herbivore��ת��ϵ��
beta=e_fh.*alpha;
mu_h=0.05;%%%%%%herbivore����Ȼ������
mu_p=0.01;%%%%%%pollinator����Ȼ������
%%
P_Nant=0.5;%%%��ʳ����Ƕ��
P_Nmut=0.5;%%%��������Ƕ��
P_Qant=0.5;%%%��������ģ��
P_Qmut=0.5;%%%��ʳ����ģ��
Cant=0.2;%%%��ʳ�������Ӷ�
Cmut=0.2;%%%�����������Ӷ�
%%
matrix_FH=get_matrix(Sf,Sh,P_Nant,P_Qant,Cant);%%%��ʳ��ϵ����
matrix_FP=get_matrix(Sf,Sp,P_Nmut,P_Qmut,Cmut);%%%���ݹ�ϵ����
%%
y0=rand(Sf+Sh+Sp,1)*10;%%%�������������ֳ�ʼֵ
t0 =0:0.1:5000;
%%
g=0.5;%%�������Ƶ�ǿ��
[t,y]=ode45(@(t,y) FHP(t,y,Sf,Sh,Sp,g,h_fh,h_fp,h_pf,r,delta_f,delta_p,alpha,beta,b_pf,c_fp,matrix_FH,matrix_FP,mu_h,mu_p),t0,y0);
figure(1)
plot(t,y)