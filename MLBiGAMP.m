function [MSE_x, MSE_h, qh]=MLBiGAMP(Input,obj)

%% load parameters
IterNum=Input.IterNum;
Kp=Input.K(1);
Kd=Input.K(2);
N1=Input.N1;
N2=Input.N2;
N3=Input.N3;
nuw2=Input.nuw2;
nuw1=Input.nuw1;
H2=obj.H2;


X_pilot=obj.X(:,1:Kp);
X_pilot_s=abs(X_pilot).^2;
Y=obj.tilde_Y;

hat_s=(H2'*H2+nuw2*eye(N2))^(-1)*H2'*Y;
hat_s1=hat_s(:,1:Kp);
hat_h1= hat_s1*X_pilot'*pinv(X_pilot*X_pilot'+N1*nuw1*eye(N1));
qh=mean(abs(hat_h1(:).^2));
hat_h_old=hat_h1;
var_vh1=ones(N2,N1)/N1;
var_vh_old=var_vh1;


%% Initialization 
Z1=zeros(N2,Kp+Kd);
Z2=zeros(N3,Kp+Kd);
V1=ones(N2,Kp+Kd);
V2=ones(N3,Kp+Kd);
V2_old=V2;
V1_old=V1;

tau1_old=ones(N2,Kp+Kd);
g1_old=zeros(N2,Kp+Kd);
tau2_old=ones(N3,Kp+Kd);
g2_old=zeros(N3,Kp+Kd);
zvar_max=1-1e-5;

hat_xData=zeros(N1,Kd);
hat_xData_old=hat_xData;
var_vxData=ones(N1,Kd);
var_vxData_old=var_vxData;
hat_s=zeros(N2,Kp+Kd);
hat_s_old=hat_s;


sqr_h2=abs(H2).^2;

MSEx_old=1;
MSEh_old=1;

MSE_x=[];
MSE_h=[];


Sigma_xData_old=ones(N1,Kd);
Sigma_h1_old=ones(N2,N1)/N1;
mes=Input.mes;
mes1=1;


%% Iteration
for kk=1:IterNum

    %% back passing
    %--second----------------------------------------- 
    [tilde_z2,tilde_vz2]=Estimator_z2(Input,obj,Z2,V2);  
    g2=(tilde_z2-Z2)./V2;
    tau2=(1./V2).*(1-min(tilde_vz2./V2,zvar_max));
    [g2,g2_old]=damping(g2,g2_old,mes);
    [tau2,tau2_old]=damping(tau2,tau2_old,mes);   

    
    Sigma_x2=(sqr_h2'*tau2).^(-1);
    R_x2=hat_s+Sigma_x2.*(H2'*g2);
    
    %--first------------------------------------------
    [tilde_z1, tilde_vz1]=Estimator_z1(Input,Z1,V1,R_x2,Sigma_x2);
    
    g1=(tilde_z1-Z1)./V1;
    tau1=(1./V1).*(1-min(tilde_vz1./V1,zvar_max));
   % [g1,g1_old]=damping(g1,g1_old,mes);
    [tau1,tau1_old]=damping(tau1,tau1_old,mes);
    

    Sigma_xData=((abs(hat_h1).^2)'*tau1(:,Kp+1:end)).^(-1);
    Gain_x1=1-Sigma_xData.*(var_vh1'*tau1(:,end-Kd+1:end));
    Gain_x1=min(1,max(0,Gain_x1));
    R_xData=hat_xData.*Gain_x1+Sigma_xData.*(hat_h1'*g1(:,end-Kd+1:end));
    [Sigma_xData,Sigma_xData_old]=damping(Sigma_xData,Sigma_xData_old,mes1);
    
    Sigma_h1=(tau1*[X_pilot_s abs(hat_xData).^2]').^(-1);
    Gain_h1=1-Sigma_h1.*(tau1*[zeros(N1,Kp) var_vxData]');
    Gain_h1=min(1,max(0,Gain_h1));
    R_h1=hat_h1.*Gain_h1+Sigma_h1.*(g1*[X_pilot hat_xData]');
    [Sigma_h1,Sigma_h1_old]=damping(Sigma_h1,Sigma_h1_old,mes1);
    
    %% Forward passing 
    %--First layer
    [hat_xData, var_vxData]=Estimator_x(obj, R_xData, Sigma_xData);
    [hat_h1, var_vh1]=Estimator_h(Input, R_h1, Sigma_h1);
    [hat_xData, hat_xData_old]=damping(hat_xData, hat_xData_old, mes);
    [hat_h1, hat_h_old]=damping(hat_h1, hat_h_old, mes);


    
     Error_X=hat_xData-obj.X_Data;
     Error_H=hat_h1-obj.H1;
     
     MSEx=norm(Error_X(:))^2/norm(obj.X_Data(:))^2;
     MSEh=norm(Error_H(:))^2/norm(obj.H1(:))^2;
    
    stop=JudgeNan(hat_xData)+JudgeNan(var_vxData)...
         +JudgeNan(hat_h1)+JudgeNan(var_vh1);
     
     if stop 
         MSE_x(kk:IterNum,1)=MSEx_old;
         MSE_h(kk:IterNum,1)=MSEh_old;
         break;
     end
    
    
     MSE_x(kk,1)=MSEx;
     MSE_h(kk,1)=MSEh;
     MSEx_old=MSEx;
     MSEh_old=MSEh;
     
     sqr_hatxData=abs(hat_xData).^2;
     sqr_hath=abs(hat_h1).^2;
     
     V1_mean=[zeros(N2,Kp) sqr_hath*var_vxData]+var_vh1*[X_pilot_s,sqr_hatxData];
     Z1_mean=hat_h1*[X_pilot hat_xData];
     V1=V1_mean+[zeros(N2,Kp) var_vh1*var_vxData];
     [V1,V1_old]=damping(V1,V1_old,mes);
%     [V1_mean,V1_mean_old]=damping(V1_mean,V1_mean_old,mes);
     Z1=Z1_mean-g1.*V1_mean;
     
     % Second layer
     [hat_s, var_vs]=Estimator_s(Input,R_x2,Sigma_x2,Z1,V1);
     [hat_s,hat_s_old]=damping(hat_s,hat_s_old,mes);
     
     
     V2=sqr_h2*var_vs;
     Z2=H2*hat_s-g2.*V2;
     [V2,V2_old]=damping(V2,V2_old,mes);
        
end
end 