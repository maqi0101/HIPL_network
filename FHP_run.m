clc;
clear; 
Sf=24;%%%%%%flower的物种数
Sh=24;%%%%%%%herbivore物种数
Sp=24;%%%%%%%%pollinator物种数
h_fh=10;%%%%%%%herbivore对flower捕食的半饱和常数
h_fp=10;%%%%%%pollinator对flower的半饱和常数
h_pf=10;%%%%%%flower对pollinator的半饱和常数
r=0.05*ones(Sf,1);%%%%%%%flower的增长率
delta_f=0.025*ones(Sf,1);%%%%%flower的密度依赖, Sf*1 matrix
delta_p=0.025*ones(Sp,1);%%%%%pollinator的密度依赖, Sf*1 matrix
alpha=unifrnd(0.2,0.3,Sf,Sh);%%%%%%(i,j)herbivore_j对flower_i的捕食系数Sf*Sh matrix
b_pf=unifrnd(0.2,0.3,Sp,Sf);%%%%%%flower从授粉获益系数
c_fp=unifrnd(0.2,0.3,Sf,Sp);%%%%%%pollinator从授粉获益系数
e_fh=0.6;%%%%%flower到herbivore的转化系数
beta=e_fh.*alpha;
mu_h=0.05;%%%%%%herbivore的自然死亡率
mu_p=0.01;%%%%%%pollinator的自然死亡率
%%
P_Nant=0.5;%%%捕食网络嵌套
P_Nmut=0.5;%%%互惠网络嵌套
P_Qant=0.5;%%%互惠网络模块
P_Qmut=0.5;%%%捕食网络模块
Cant=0.2;%%%捕食网络连接度
Cmut=0.2;%%%互惠网络连接度
%%
matrix_FH=get_matrix(Sf,Sh,P_Nant,P_Qant,Cant);%%%捕食关系矩阵
matrix_FP=get_matrix(Sf,Sp,P_Nmut,P_Qmut,Cmut);%%%互惠关系矩阵
%%
y0=rand(Sf+Sh+Sp,1)*10;%%%列向量，各物种初始值
t0 =0:0.1:5000;
%%
g=0.5;%%传粉限制的强度
[t,y]=ode45(@(t,y) FHP(t,y,Sf,Sh,Sp,g,h_fh,h_fp,h_pf,r,delta_f,delta_p,alpha,beta,b_pf,c_fp,matrix_FH,matrix_FP,mu_h,mu_p),t0,y0);
figure(1)
plot(t,y)