function [hatx_old, MSEx]=SecondPhase(Input, obj, hath)
Kd=Input.K(2);
P=Input.P;
M=Input.M;
N=Input.N;
C=obj.C;
IterNum=Input.IterNum;
nuw1=Input.nuw1;
mes=Input.mes;



%% Parameters initialization 
Z2=zeros(P,Kd);
V2=ones(P,Kd);
Z1=zeros(M,Kd);
V1=ones(M,Kd);
V1_old=V1;
Z1_old=Z1;
V2_old=V2;
Z2_old=Z2;
g2_old=zeros(P, Kd);
tau2_old=ones(P, Kd);
g1_old=zeros(M, Kd);
tau1_old=ones(M, Kd);


hats=zeros(M,Kd);
hatx=zeros(N,Kd);
hatx_old=hatx;

sqrH=abs(hath).^2;
sqrtH=sqrH';
Ht=hath';


sqrC=abs(C).^2;
sqrtC=sqrC';
Ct=C';


MSE_x_old=1;
for kk=1:IterNum
    [hatz2, varz2]=EstZ_Second(Input, obj, Z2, V2);
    g2=(hatz2-Z2)./V2;
    tau2=(V2-varz2)./(V2.^2);
    [g2, g2_old]=damping(g2, g2_old, mes);
    [tau2, tau2_old]=damping(tau2, tau2_old, mes);
    
    Sigmas=(sqrtC*tau2).^(-1);
    Rs=hats+Sigmas.*(Ct*g2);
    
    varz1=1./(1./(nuw1+Sigmas)+1./V1);
    hatz1=varz1.*(Rs./(nuw1+Sigmas)+Z1./V1);
    
    g1=(hatz1-Z1)./V1;
    tau1=(V1-varz1)./(V1.^2);
    [g1, g1_old]=damping(g1, g1_old, mes);
    [tau1, tau1_old]=damping(tau1, tau1_old, mes);
    
    Sigmax=(sqrtH*tau1).^(-1);
    Rx=hatx+Sigmax.*(Ht*g1);
    
    [hatx, varx]=EstimatorX(obj, Rx, Sigmax);
    
    error_X=hatx-obj.X_Data;
    MSE_x=norm(error_X(:))^2/(norm(obj.X_Data(:))^2);
    
    flag=isnan(MSE_x);
    
    if flag
        MSEx(kk:IterNum,1)=MSE_x_old;
        break;
    end 
    hatx_old=hatx;
    MSEx(kk,1)=MSE_x;
    MSE_x_old=MSE_x;
    
    V1=sqrH*varx;
    Z1=hath*hatx-g1.*V1;
    [V1, V1_old]=damping(V1, V1_old, mes);
    [Z1, Z1_old]=damping(Z1, Z1_old, mes);
    
    vars=1./(1./Sigmas+1./(nuw1+V1));
    hats=vars.*(Rs./Sigmas+Z1./(nuw1+V1));
    
    V2=sqrC*vars;
    Z2=C*hats-g2.*V2;
    
    [V2, V2_old]=damping(V2, V2_old, mes);
    [Z2, Z2_old]=damping(Z2, Z2_old, mes);
end