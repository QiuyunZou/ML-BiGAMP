function flag=MLBiGAMP_MC(Input,obj)
        

%% load parameters
N1=Input.N1;
N2=Input.N2;
N3=Input.N3;
K=Input.K;
IterNum=Input.IterNum;
omega=obj.omega;
step=0.5;
stepMin=0.05;
stepMax=0.65;
stepIncr=1.1;
step1=1;

zvar_max=0.99;

hath=randn(N2, N1);
var_vh=10*ones(N2, N1);

hatx1=randn(N1, K);
var_vx1=10*ones(N1, K);

hath_old=zeros(N2, N1);
hatx_old=zeros(N1, K);

g1=zeros(N2,K);
g1_old=g1;
tau1_old=zeros(N2,K);

V1_old=0;
V1_mean_old=0;

H2=obj.H2;
sqr_h2=abs(H2).^2;

g2_old=zeros(N3, K);
tau2_old=ones(N3,K);
hat_x2=randn(N2,K);
hat_x2_old=randn(N2,K);

Z2=zeros(N3, K);
V2=ones(N3, K);
V2_old=V2;

Z1=zeros(N2, K);
V1=ones(N2,K);

MSEz=zeros(IterNum, 1);


maskMatrix=zeros(N2,K);
maskMatrix(omega)=1;

flag=false;
for kk=1:IterNum
    

    
%% Back passing, second layer
[tilde_z2,tilde_vz2]=Estimator_z2(Input,obj,Z2,V2); 
g2=(tilde_z2-Z2)./V2;
tau2=(1-min(tilde_vz2./V2,zvar_max))./V2;
[g2,g2_old]=damping(g2,g2_old, step1);
[tau2,tau2_old]=damping(tau2,tau2_old, step1);   

Sigma_x2=(sqr_h2'*tau2).^(-1);
R_x2=hat_x2+Sigma_x2.*(H2'*g2);
    
[tilde_z1, tilde_vz1]=Estimator_z1(Input,Z1,V1,R_x2,Sigma_x2);
tilde_vz1=real(tilde_vz1);

VInvMask=(1./V1).* maskMatrix;
g1 = VInvMask.*(tilde_z1-Z1);
tau1 = VInvMask.*(1-min(tilde_vz1./V1,zvar_max));
[g1,g1_old]=damping(g1,g1_old,step1);
[tau1,tau1_old]=damping(tau1,tau1_old,step1);


if kk>1 
    step1=step;
end

step = max(stepMin, stepIncr*step);   % increase damping
step=  min(step, stepMax);   
    

Sigma_x=1./((abs(hath).^2)'*tau1+eps);
Gain_x =1-Sigma_x.*(var_vh'*tau1);
Gain_x =min(1,max(0,Gain_x));
R_x=hatx1.*Gain_x+Sigma_x.*(hath'*g1);

Sigma_h=(tau1*(abs(hatx1).^2)').^(-1);
Gain_h =1-Sigma_h.*(tau1*var_vx1');
Gain_h = min(1,max(0,Gain_h));
R_h=hath.*Gain_h+Sigma_h.*(g1*hatx1');

%% ----
[hatx1,var_vx1]=Estimator_x(R_x, Sigma_x);
[hath,var_vh]=Estimator_h(R_h,Sigma_h);
[hatx1,hatx_old]  = damping(hatx1,hatx_old,step1);
[hath,hath_old]  = damping(hath,hath_old,step1);


sqr_hatx=abs(hatx1).^2;
sqr_hath=abs(hath).^2;

V1_mean=var_vh*sqr_hatx+ sqr_hath*var_vx1;
Z1_mean=hath*hatx1;
V1=V1_mean+ var_vh*var_vx1;
[V1_mean, V1_mean_old]=damping(V1_mean, V1_mean_old, step1);
[V1, V1_old]=damping(V1, V1_old, step1);
Z1=Z1_mean-g1.*V1_mean;


mse=10*log10(norm(Z1_mean - obj.Z1,'fro')^2 / norm(obj.Z1,'fro')^2);
MSEz(kk)=mse;

if mse<-50 
    flag=true;
    break;
elseif JudgeNan(Z1_mean)
    break;
end

[hat_x2, var_vx2]  =Estimator_s(Input,obj, R_x2,Sigma_x2,Z1,V1);
[hat_x2,hat_x2_old]=damping(hat_x2,hat_x2_old, step1);


V2=sqr_h2*var_vx2;
Z2=H2*hat_x2-g2.*V2;
[V2,V2_old]=damping(V2,V2_old,step1);


end


end
