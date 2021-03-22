function MSE_error=MLVAMP(Input, obj)

N1=Input.N1;
N2=Input.N2;
N3=Input.N3;
N4=Input.N4;

H1=obj.H1;
H2=obj.H2;
H3=obj.H3;
mes=Input.mes;
IterNum=Input.IterNum;

u1_plus=ones(N2,1);
u1_plus_inv=1./u1_plus;
m1_plus=zeros(N2,1);
u1_plus_inv_old=u1_plus_inv;
m1_plus_old=m1_plus;

u2_plus=ones(N3,1);
u2_plus_inv=1./u2_plus;
m2_plus=zeros(N3,1);
u2_plus_inv_old=u2_plus_inv;
m2_plus_old=m2_plus;


u3_plus=ones(N4,1);
u3_plus_inv=1./u3_plus;
m3_plus=zeros(N4,1);
u3_plus_inv_old=u3_plus_inv;
m3_plus_old=m3_plus;

Sigma1_plus=ones(N1,1);
Sigma1_plus_inv=1./Sigma1_plus;
R1_plus=zeros(N1,1);
Sigma1_plus_inv_old=Sigma1_plus_inv;
R1_plus_old=R1_plus;

Sigma2_plus=ones(N2,1);
Sigma2_plus_inv=1./Sigma2_plus;
R2_plus=zeros(N2,1);
Sigma2_plus_inv_old=Sigma2_plus_inv;
R2_plus_old=R2_plus;


Sigma3_plus=ones(N3,1);
Sigma3_plus_inv=1./Sigma3_plus;
R3_plus=zeros(N3,1);
Sigma3_plus_inv_old=Sigma3_plus_inv;
R3_plus_old=R3_plus;

MSE_error=zeros(IterNum,1);
MSE_old=1;

for ii=1:IterNum
    
    % 3th layer
    [hatz3_sub,vz3_sub]=Estimator_z3(Input, obj, m3_plus./u3_plus_inv,1./u3_plus_inv);
    u3_sub=max(vz3_sub./(1-vz3_sub.*u3_plus_inv),eps);          
    m3_sub=u3_sub.*(hatz3_sub./vz3_sub-m3_plus);
        
    Qx3_sub=(H3'*diag(1./u3_sub)*H3+diag(Sigma3_plus_inv))^(-1);
    hatx3_sub=Qx3_sub*(H3'*diag(1./u3_sub)*m3_sub+R3_plus);                     
    vx3_sub=real(diag(Qx3_sub));
    Sigma3_sub=max(vx3_sub./(1-vx3_sub.*Sigma3_plus_inv),eps);
    R3_sub=Sigma3_sub.*(hatx3_sub./vx3_sub-R3_plus);
    
    % 2ed layer
    [hatz2_sub,vz2_sub]=Estimator_z2(Input, m2_plus./u2_plus_inv,1./u2_plus_inv, R3_sub, Sigma3_sub);
    u2_sub=max(vz2_sub./(1-vz2_sub.*u2_plus_inv),eps);          
    m2_sub=u2_sub.*(hatz2_sub./vz2_sub-m2_plus);
        
    Qx2_sub=(H2'*diag(1./u2_sub)*H2+diag(Sigma2_plus_inv))^(-1);
    hatx2_sub=Qx2_sub*(H2'*diag(1./u2_sub)*m2_sub+R2_plus);                     
    vx2_sub=real(diag(Qx2_sub));
    Sigma2_sub=max(vx2_sub./(1-vx2_sub.*Sigma2_plus_inv),eps);
    R2_sub=Sigma2_sub.*(hatx2_sub./vx2_sub-R2_plus);
    
    % 1st layer
    [hatz1_sub,vz1_sub]=Estimator_z1(Input, m1_plus./u1_plus_inv,1./u1_plus_inv, R2_sub, Sigma2_sub);
    u1_sub=max(vz1_sub./(1-vz1_sub.*u1_plus_inv),eps);         
    m1_sub=u1_sub.*(hatz1_sub./vz1_sub-m1_plus);
        
    Qx1_sub=(H1'*diag(1./u1_sub)*H1+diag(Sigma1_plus_inv))^(-1);
    hatx1_sub=Qx1_sub*(H1'*diag(1./u1_sub)*m1_sub+R1_plus);                    
    vx1_sub=real(diag(Qx1_sub));
    Sigma1_sub=max(vx1_sub./(1-vx1_sub.*Sigma1_plus_inv),eps);
    R1_sub=Sigma1_sub.*(hatx1_sub./vx1_sub-R1_plus);
    
    %% Forward passing
    [hatx1_plus,vx1_plus]=Estimator_x1(obj, R1_sub, Sigma1_sub);
    MSE=mean(abs(hatx1_plus(:)-obj.X(:)).^2);
    MSE_error(ii,1)=MSE;
    
    
    if isnan(MSE) || isinf(MSE) 
        MSE_error(ii:end,1)=MSE_old;
        break;
    end
    MSE_old=MSE;
    
    vx1_plus=max(vx1_plus,5e-13);
    Sigma1_plus_inv=(Sigma1_sub-vx1_plus)./vx1_plus./Sigma1_sub;
    R1_plus=(hatx1_plus.*Sigma1_sub-R1_sub.*vx1_plus)./vx1_plus./Sigma1_sub;
    
    negldx=Sigma1_plus_inv<eps;
    Sigma1_plus_inv(negldx)=Sigma1_plus_inv_old(negldx);
    R1_plus(negldx)=R1_plus_old(negldx);
    
    % damping
    [Sigma1_plus_inv,Sigma1_plus_inv_old]=damping(Sigma1_plus_inv,Sigma1_plus_inv_old, mes);
    [R1_plus,R1_plus_old]=damping(R1_plus,R1_plus_old, mes);
      
    Qx1_plus=(H1'*diag(1./u1_sub)*H1+diag(Sigma1_plus_inv))^(-1);
    mx1_plus=Qx1_plus*(H1'*diag(1./u1_sub)*m1_sub+R1_plus);
    hatz1_plus=H1*mx1_plus;
    vz1_plus=real(diag(H1*Qx1_plus*H1'));
          
    vz1_plus=max(vz1_plus,5e-13);
    u1_plus_inv=(u1_sub-vz1_plus)./u1_sub./vz1_plus;
    m1_plus=(hatz1_plus.*u1_sub-m1_sub.*vz1_plus)./u1_sub./vz1_plus;
    
    %% �޳���Ԫ�� 
    negldx=u1_plus_inv<0;
    u1_plus_inv(negldx)=u1_plus_inv_old(negldx);
    m1_plus(negldx)=m1_plus_old(negldx);
     
     %% damping
    [u1_plus_inv,u1_plus_inv_old]=damping(u1_plus_inv,u1_plus_inv_old,mes);
    [m1_plus,m1_plus_old]=damping(m1_plus,m1_plus_old,mes);
     
    
    % 2ed layer
    [hatx2_plus,vx2_plus]=Estimator_x2(Input, m1_plus./u1_plus_inv,1./u1_plus_inv, R2_sub, Sigma2_sub);
    vx2_plus=max(vx2_plus,5e-13);
    Sigma2_plus_inv=(Sigma2_sub-vx2_plus)./vx2_plus./Sigma2_sub;
    R2_plus=(hatx2_plus.*Sigma2_sub-R2_sub.*vx2_plus)./vx2_plus./Sigma2_sub;
    
    negldx=Sigma2_plus_inv<eps;
    Sigma2_plus_inv(negldx)=Sigma2_plus_inv_old(negldx);
    R2_plus(negldx)=R2_plus_old(negldx);
    
    % damping
    [Sigma2_plus_inv,Sigma2_plus_inv_old]=damping(Sigma2_plus_inv,Sigma2_plus_inv_old, mes);
    [R2_plus,R2_plus_old]=damping(R2_plus,R2_plus_old, mes);
      
    Qx2_plus=(H2'*diag(1./u2_sub)*H2+diag(Sigma2_plus_inv))^(-1);
    mx2_plus=Qx2_plus*(H2'*diag(1./u2_sub)*m2_sub+R2_plus);
    hatz2_plus=H2*mx2_plus;
    vz2_plus=real(diag(H2*Qx2_plus*H2'));
          
    vz2_plus=max(vz2_plus,5e-13);
    u2_plus_inv=(u2_sub-vz2_plus)./u2_sub./vz2_plus;
    m2_plus=(hatz2_plus.*u2_sub-m2_sub.*vz2_plus)./u2_sub./vz2_plus;
    
    negldx=u2_plus_inv<0;
    u2_plus_inv(negldx)=u2_plus_inv_old(negldx);
    m2_plus(negldx)=m2_plus_old(negldx);
     
     %% damping
    [u2_plus_inv,u2_plus_inv_old]=damping(u2_plus_inv,u2_plus_inv_old,mes);
    [m2_plus,m2_plus_old]=damping(m2_plus,m2_plus_old,mes);
    
    % 3th layer
    [hatx3_plus,vx3_plus]=Estimator_x3(Input, m2_plus./u2_plus_inv,1./u2_plus_inv, R3_sub, Sigma3_sub);
    vx3_plus=max(vx3_plus,5e-13);
    Sigma3_plus_inv=(Sigma3_sub-vx3_plus)./vx3_plus./Sigma3_sub;
    R3_plus=(hatx3_plus.*Sigma3_sub-R3_sub.*vx3_plus)./vx3_plus./Sigma3_sub;
    
    negldx=Sigma3_plus_inv<eps;
    Sigma3_plus_inv(negldx)=Sigma3_plus_inv_old(negldx);
    R3_plus(negldx)=R3_plus_old(negldx);
    
    % damping
    [Sigma3_plus_inv,Sigma3_plus_inv_old]=damping(Sigma3_plus_inv,Sigma3_plus_inv_old, mes);
    [R3_plus,R3_plus_old]=damping(R3_plus,R3_plus_old, mes);
      
    Qx3_plus=(H3'*diag(1./u3_sub)*H3+diag(Sigma3_plus_inv))^(-1);
    mx3_plus=Qx3_plus*(H3'*diag(1./u3_sub)*m3_sub+R3_plus);
    hatz3_plus=H3*mx3_plus;
    vz3_plus=real(diag(H3*Qx3_plus*H3'));
          
    vz3_plus=max(vz3_plus,5e-13);
    u3_plus_inv=(u3_sub-vz3_plus)./u3_sub./vz3_plus;
    m3_plus=(hatz3_plus.*u3_sub-m3_sub.*vz3_plus)./u3_sub./vz3_plus;
    
    negldx=u3_plus_inv<0;
    u3_plus_inv(negldx)=u3_plus_inv_old(negldx);
    m3_plus(negldx)=m3_plus_old(negldx);
     
     %% damping
    [u3_plus_inv,u3_plus_inv_old]=damping(u3_plus_inv,u3_plus_inv_old,mes);
    [m3_plus,m3_plus_old]=damping(m3_plus,m3_plus_old,mes);
end
end

