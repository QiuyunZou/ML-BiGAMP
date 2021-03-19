function [ hath_old, MSEh]=FirstPhase(Input, obj)

%% Load parameters
Kp=Input.K(1);
P=Input.P;
M=Input.M;
N=Input.N;
Xp=obj.X_pilot;
C=obj.C;
IterNum=Input.IterNum;
nuw1=Input.nuw1;
mes=0.8;



%% Initialize parameters
Z2=zeros(P,Kp);
V2=ones(P,Kp);
Z1=zeros(M,Kp);
V1=ones(M,Kp);
V1_old=V1;
Z1_old=Z1;
V2_old=V2;
Z2_old=Z2;

sqrC=abs(C).^2;
sqrtC=sqrC';
Ct=C';


hats=zeros(M,Kp);
hath=zeros(M,N);

sqr_X=abs(Xp).^2;

MSE_h_old=1;
flag=0;
%% Load parameters
for ii=1:IterNum
    [hatz2, varz2]=EstZ_First(Input, obj, Z2, V2);
    g2=(hatz2-Z2)./V2;
    tau2=(V2-varz2)./(V2.^2);
    
    Sigmas=(sqrtC*tau2).^(-1);
    Rs=hats+Sigmas.*(Ct*g2);
    
    varz1=1./(1./(nuw1+Sigmas)+1./V1);
    hatz1=varz1.*(Rs./(nuw1+Sigmas)+Z1./V1);
    
    g1=(hatz1-Z1)./V1;
    tau1=(V1-varz1)./(V1.^2);
    
    Sigma_h=(tau1*(abs(Xp).^2)').^(-1);
    R_h=hath+Sigma_h.*(g1*Xp');
    
    [hath, varh]=Estimator_h(Input, R_h, Sigma_h);
    
    
    
    Error_H=hath-obj.H;
    MSE_h=norm(Error_H(:))^2/(norm(obj.H(:))^2);
    
    if JudgeNan(hath) ||  JudgeNan(varh)
        MSEh(ii:IterNum,1)=MSE_h_old;
        break;
    end
    
    hath_old=hath;
    MSE_h_old=MSE_h;
    MSEh(ii,1)=MSE_h;
    
    V1=varh*sqr_X;
    Z1=hath*Xp-g1.*V1;
    
    [V1, V1_old]=damping(V1, V1_old, mes);
    [Z1, Z1_old]=damping(Z1, Z1_old, mes);
    
    vars=1./(1./Sigmas+1./(nuw1+V1));
    hats=vars.*(Rs./Sigmas+Z1./(nuw1+V1));
    
    V2=sqrC*vars;
    Z2=C*hats-g2.*V2;
    
    [V2, V2_old]=damping(V2, V2_old, mes);
    [Z2, Z2_old]=damping(Z2, Z2_old, mes);
    
    
end

end 