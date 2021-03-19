function [hat_xData, MSE_x, MSE_h]=JCD(Input,obj)

%% load parameters
IterNum=Input.IterNum;
K1=Input.K(1);
K2=Input.K(2);
N=Input.N;
M=Input.M;
P=Input.P;
nuw2=Input.nuw2;
nuw1=Input.nuw1;
C=obj.C;


X_pilot=obj.X(:,1:K1);
X_pilot_s=abs(X_pilot).^2;
Y=obj.tilde_Y;

hat_s=(C'*C+nuw2*eye(M))^(-1)*C'*Y;
hat_s1=hat_s(:,1:K1);
hat_h= hat_s1*X_pilot'*pinv(X_pilot*X_pilot'+N*nuw1*eye(N));
hat_h_old=hat_h;
var_vh=ones(M,N)/N;


%% Initialization 
Z1=zeros(M,K1+K2);
Z2=zeros(P,K1+K2);
V1=ones(M,K1+K2);
V2=ones(P,K1+K2);
V2_old=V2;
V1_old=V1;

tau1_old=ones(M,K1+K2);
tau2_old=ones(P,K1+K2);
g2_old=zeros(P,K1+K2);
zvar_max=1-1e-5;

hat_xData=zeros(N,K2);
hat_xData_old=hat_xData;
var_vxData=ones(N,K2);
hat_s=zeros(M,K1+K2);
g1_old=zeros(M,K1+K2);
hat_s_old=hat_s;


sqr_c=abs(C).^2;


Sigma_xData_old=ones(N,K2);
Sigma_h_old=ones(M,N)/N;
Z2_old=zeros(P,K1+K2);
Z1_old=Z1;
mes=Input.mes;
mes1=1;


MSEx_old=1;
MSEh_old=1;
%% Iteration
for kk=1:IterNum

    %% back passing
    %--second----------------------------------------- 
    [tilde_z2,tilde_vz2]=EstZ_Second(Input, obj, Z2, V2);  
    g2=(tilde_z2-Z2)./V2;
    tau2=(1./V2).*(1-min(tilde_vz2./V2,zvar_max));
    [g2,g2_old]=damping(g2,g2_old,mes);
    [tau2,tau2_old]=damping(tau2,tau2_old,mes);   

    
    Sigma_s=(sqr_c'*tau2).^(-1);
    R_s=hat_s+Sigma_s.*(C'*g2);
    
    %--first------------------------------------------
    [tilde_z1, tilde_vz1]=Estimator_z1(Input,Z1,V1,R_s,Sigma_s);
    
    g1=(tilde_z1-Z1)./V1;
    tau1=(1./V1).*(1-min(tilde_vz1./V1,zvar_max));
    [g1, g1_old]=damping(g1, g1_old, mes);
    [tau1,tau1_old]=damping(tau1,tau1_old,mes);
    

    Sigma_xData=((abs(hat_h).^2)'*tau1(:,K1+1:end)).^(-1);
    Gain_x=1-Sigma_xData.*(var_vh'*tau1(:,end-K2+1:end));
    Gain_x=min(1,max(0,Gain_x));
    R_xData=hat_xData.*Gain_x+Sigma_xData.*(hat_h'*g1(:,end-K2+1:end));
    [Sigma_xData,Sigma_xData_old]=damping(Sigma_xData,Sigma_xData_old,mes1);
    
    Sigma_h=(tau1*[X_pilot_s abs(hat_xData).^2]').^(-1);
    Gain_h=1-Sigma_h.*(tau1*[zeros(N,K1) var_vxData]');
    Gain_h=min(1,max(0,Gain_h));
    R_h=hat_h.*Gain_h+Sigma_h.*(g1*[X_pilot hat_xData]');
    [Sigma_h,Sigma_h_old]=damping(Sigma_h,Sigma_h_old,mes1);
    
    %% Forward passing 
    %--First layer
    [hat_xData, var_vxData]=Estimator_x(obj,R_xData,Sigma_xData);
    [hat_h,var_vh]=Estimator_h(Input,R_h,Sigma_h);
    

    
    MSEx=norm(hat_xData(:)-obj.X_Data(:))^2/(norm(obj.X_Data(:))^2);
    MSE_x(kk,1)=MSEx;

    H=obj.H; 
    MSEh=norm(hat_h(:)-H(:))^2/(norm(H(:))^2);
    MSE_h(kk,1)=MSEh;
     
     
     if JudgeNan(hat_xData) || JudgeNan(var_vxData)
         hat_xData=hat_xData_old;
         MSE_x(kk:IterNum,1)=MSEx_old;
         MSE_h(kk:IterNum,1)=MSEh_old;
         break;
     end

    
     MSEh_old=MSEh;
     MSEx_old=MSEx;
     
     
     sqr_hatxData=abs(hat_xData).^2;
     sqr_hath=abs(hat_h).^2;
     
     [hat_xData,hat_xData_old]=damping(hat_xData,hat_xData_old,mes);
     [hat_h,hat_h_old]=damping(hat_h,hat_h_old,mes);
     
     V1_mean=[zeros(M,K1) sqr_hath*var_vxData]+var_vh*[X_pilot_s,sqr_hatxData];
     Z1_mean=hat_h*[X_pilot hat_xData];
     V1=V1_mean+[zeros(M,K1) var_vh*var_vxData];
     [V1,V1_old]=damping(V1,V1_old,mes);
     Z1=Z1_mean-g1.*V1_mean;
%     [Z1, Z1_old]=damping(Z1, Z1_old, mes);
     
     % Second layer
     [hat_s, var_vs]=Estimator_s(Input,R_s,Sigma_s,Z1,V1);
   %  [hat_s,hat_s_old]=damping(hat_s,hat_s_old,mes);
     
     
     V2=sqr_c*var_vs;
     Z2=C*hat_s-g2.*V2;
     [V2,V2_old]=damping(V2,V2_old,mes);
  %   [Z2, Z2_old]=damping(Z2, Z2_old, mes);
        
end
end 