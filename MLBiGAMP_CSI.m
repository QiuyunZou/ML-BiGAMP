function MSE_array=MLBiGAMP_CSI(Input,obj)

%% load parameters
IterNum=Input.IterNum;
N1=Input.N1;
N2=Input.N2;
N3=Input.N3;
N4=Input.N4;

H1=obj.H1;
H2=obj.H2;
H3=obj.H3;


%% Initialization 
Z1=zeros(N2, 1);
V1=ones(N2, 1);
V1_old=V1;

Z2=zeros(N3, 1);
V2=ones(N3, 1);
V2_old=V2;

Z3=zeros(N4, 1);
V3=ones(N4, 1);
V3_old=V3;

zvar_max=1-1e-5;

hat_x1=zeros(N1, 1);
hat_x2=zeros(N2, 1);
hat_x3=zeros(N3, 1);


sqr_H1=abs(H1).^2;
sqr_H2=abs(H2).^2;
sqr_H3=abs(H3).^2;

mes=Input.mes;
MSEx_old=1;


%% Iteration
for kk=1:IterNum
    %% back passing
    [tilde_z3,tilde_vz3]=Estimator_z3(Input, obj, Z3, V3); 
    g3=(tilde_z3-Z3)./V3;
    tau3=(1./V3).*(1-min(tilde_vz3./V3,zvar_max));
    Sigma_x3=(sqr_H3'*tau3).^(-1);
    R_x3=hat_x3+Sigma_x3.*(H3'*g3);
    
    % 2-ed layer
    [tilde_z2,tilde_vz2]=Estimator_z2(Input, Z2, V2, R_x3, Sigma_x3); 
    g2=(tilde_z2-Z2)./V2;
    tau2=(1./V2).*(1-min(tilde_vz2./V2,zvar_max));
    Sigma_x2=(sqr_H2'*tau2).^(-1);
    R_x2=hat_x2+Sigma_x2.*(H2'*g2);

   
    % 1-st layer
    [tilde_z1,tilde_vz1]=Estimator_z1(Input, Z1, V1, R_x2, Sigma_x2); 
    g1=(tilde_z1-Z1)./V1;
    tau1=(1./V1).*(1-min(tilde_vz1./V1,zvar_max));
    Sigma_x1=(sqr_H1'*tau1).^(-1);
    R_x1=hat_x1+Sigma_x1.*(H1'*g1);

    %% forward passing 
    [hat_x1, var_vx1]=Estimator_x1(obj, R_x1, Sigma_x1);
    MSEx=mean(abs(hat_x1(:)-obj.X(:)).^2);
    MSE_array(kk,1)=MSEx;
    
    stop=JudgeNan(MSEx);  
    
     if stop>=1 || MSEx>MSEx_old
         MSE_array(kk:IterNum,1)=MSEx_old;
         break;
     end
     MSEx_old=MSEx;
     
     V1=sqr_H1*var_vx1;
     Z1=H1*hat_x1-g1.*V1;
     [V1,V1_old]=damping(V1,V1_old,mes);
     
     [hat_x2, var_vx2]=Estimator_x2(Input, Z1, V1, R_x2, Sigma_x2); 
     
     V2=sqr_H2*var_vx2;
     Z2=H2*hat_x2-g2.*V2;
     [V2,V2_old]=damping(V2,V2_old,mes);
     
     [hat_x3, var_vx3]=Estimator_x3(Input, Z2, V2, R_x3, Sigma_x3); 
     
     V3=sqr_H3*var_vx3;
     Z3=H3*hat_x3-g3.*V3;
     [V3,V3_old]=damping(V3,V3_old,mes);
        
end
end 