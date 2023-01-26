function dydt = FHP(t,y,Sf,Sh,Sp,g,h_fh,h_fp,h_pf,r,delta_f,delta_p,alpha,beta,b_pf,c_fp,matrix_FH,matrix_FP,mu_h,mu_p)
dydt=zeros(Sf+Sh+Sp,1);%Sf:flower;Sh:herbivore;Sp:pollinator，列向量
Fi=y(1:Sf);% flower
Hi=y(Sf+1:Sf+Sh);%herbivore
Pi=y(Sf+Sh+1:Sf+Sh+Sp);%pollinator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gpi=ones(Sp,1);%食草动物对植物i的干扰作用
Gpi=1./(1+g*matrix_FH*Hi);%Sp*1 matrix
%%
%食草动物和植物间相互作用
F_fh=matrix_FH.*(Fi*Hi')./(h_fh+repmat((matrix_FH'*Fi)',[Sf,1]));%check
%%
%传粉者从植物获得效应
E1_fp=matrix_FP.*((Gpi.*Fi)*Pi')./(h_fp+repmat((matrix_FP'*Fi)',[Sf,1]));%check
%%
%植物从传粉者获得效应
E2_pf=matrix_FP'.*(Pi*(Gpi.*Fi)')./(h_pf+repmat((matrix_FP*Pi)',[Sp,1]));%check
%%
dFi=Fi.*(r-delta_f.*Fi)-sum(alpha.*F_fh,2)+(sum(b_pf.*E2_pf))';%check
dHi=(sum(beta.*F_fh))'-mu_h.*Hi;%check
dPi=(sum(c_fp.*E1_fp))'-delta_p.*Pi.^2-mu_p.*Pi;%check
dydt(1:Sf)=dFi;
dydt(Sf+1:Sf+Sh)=dHi;
dydt(Sf+Sh+1:Sf+Sh+Sp)=dPi;