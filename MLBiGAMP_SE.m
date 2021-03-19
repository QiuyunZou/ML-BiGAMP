function [MSE_x,MSE_h]=MLBiGAMP_SE(Input,q_h)

N1=Input.N1;
N2=Input.N2;
N3=Input.N3;
Kp=Input.K(1);
Kd=Input.K(2);
K=Kp+Kd;
rho=Kp/K;

nuw1=Input.nuw1;
nuw2=Input.nuw2;
IterNum=Input.IterNum;



n_bit=Input.n_bit;
Gaussian=@(x,m,v) 1./sqrt(2*pi*v).*exp(-(x-m).^2./(2*v));
if Input.ADC_switch==1
    switch n_bit
        case 1
            quan_step=2;
        case 2
             quan_step=sqrt(0.25);
        case 3
            quan_step=sqrt(0.25);
        otherwise
            quan_step=min(max(abs(Y_real)),max(abs(Y_imag)))/2^(n_bit-1);    
    end
    Q_out = (-2^n_bit+1:2:2^n_bit-1)*quan_step/2 ;
    y_up=Q_out+quan_step/2;
    y_low=Q_out-quan_step/2;
    y_up(end)=Inf;
    y_low(1)=-Inf;
end

Tx1=1;
Th1=1/N1;
Th2=1/N2;
qh2=Th2;
Tz1=N1*Tx1*Th1;
Tx2=Tz1+nuw1;
Tz2=N2*Tx2*Th2;
 eta=@(y,V,xi) (y-sqrt((Tz2-V)/2).*xi)./sqrt((nuw2+V)/2);

V2_P=1;
V2_D=1;

V1_P=1;
V1_D=1;

qx1=0;
qh1=1/(1.5*N1);
%qh=mean(q_h);

tau2_P_old=0;
tau2_D_old=0;

tau1_P_old=0;
tau1_D_old=0;

V1_P_old=1;
V1_D_old=1;
V2_P_old=1;
V2_D_old=1;
Sigmax_D_old=1;
Sigmah_old=1/N1;


MSEx_old=1;
MSEh_old=1/N1;

mes1=1;
mes=Input.mes;

Sigmax_D=0.75;
Sigmah=1/(2*N1);
Sigmas_P=1;
Sigmas_D=1;
for kk=1:IterNum
%% Forward passing
    gamma=1/Sigmax_D;
    MSEx=1-integral(@(z) tanh(gamma+sqrt(gamma)*z).*Gaussian(z,0,1),-Inf,Inf);
    MSEh=1/(1/Sigmah+N1);
    %[MSEx,MSEx_old]=damping(MSEx,MSEx_old,mes1);
    %[MSEh,MSEh_old]=damping(MSEh,MSEh_old,mes1);
    
    MSE_x(kk,1)=MSEx;
    MSE_h(kk,1)=MSEh*N1;
    MSEx_old=MSEx;
    MSEh_old=MSEh;
    
    qx1=1-MSEx;
    qh1=1/N1-MSEh;   
    
    V1_P=N1*(Tx1-qh1);
    V1_D=N1*(Tx1*Th1-qx1*qh1);
    [V1_P,V1_P_old]=damping(V1_P,V1_P_old,mes);
    [V1_D,V1_D_old]=damping(V1_D,V1_D_old,mes);
    
    vs_P=Sigmas_P*(nuw1+V1_P)/(Sigmas_P+(nuw1+V1_P));
    vs_D=Sigmas_D*(nuw1+V1_D)/(Sigmas_D+(nuw1+V1_D));
    
    V2_P=N2*Th2*vs_P;
    V2_D=N2*Th2*vs_D;
    [V2_P,V2_P_old]=damping(V2_P,V2_P_old,mes);
    [V2_D,V2_D_old]=damping(V2_D,V2_D_old,mes);
    
    %% Back passing
     if Input.ADC_switch==1  %ADC case
        TemP=0;
        TemD=0;
        if V2_P>Tz2
            V2_P=Tz2-eps;
        end
        for ii=1:length(Q_out)
            TemP=TemP+integral(@(xi) (normpdf(eta(y_up(ii),V2_P,xi))-normpdf(eta(y_low(ii),V2_P,xi))).^2 ...
            ./(normcdf(eta(y_up(ii),V2_P,xi))-normcdf(eta(y_low(ii),V2_P,xi))+eps).*Gaussian(xi,0,1),-Inf,Inf);
        end
    
        for ii=1:length(Q_out)
            TemD=TemD+integral(@(xi) (normpdf(eta(y_up(ii),V2_D,xi))-normpdf(eta(y_low(ii),V2_D,xi))).^2 ...
            ./(normcdf(eta(y_up(ii),V2_D,xi))-normcdf(eta(y_low(ii),V2_D,xi))+eps).*Gaussian(xi,0,1),-Inf,Inf);
        end
        tildevz2_P=V2_P-V2_P^2/(V2_P+nuw2)*TemP;
        tildevz2_D=V2_D-V2_D^2/(V2_D+nuw2)*TemD;
     else
        tildevz2_P=V2_P*nuw2/(nuw2+V2_P);
        tildevz2_D=V2_D*nuw2/(nuw2+V2_D);
     end
    
    tau2_P=1/V2_P*(1-tildevz2_P/V2_P);
    tau2_D=1/V2_D*(1-tildevz2_D/V2_D);
    [tau2_P,tau2_P_old]=damping(tau2_P,tau2_P_old,mes);
    [tau2_D,tau2_D_old]=damping(tau2_D,tau2_D_old,mes);
    
    Sigmas_P=(N3*qh2*tau2_P)^(-1);
    Sigmas_D=(N3*qh2*tau2_D)^(-1);
    
    tildevz1_P=V1_P*(Sigmas_P+nuw1)/(V1_P+(Sigmas_P+nuw1));
    tildevz1_D=V1_D*(Sigmas_D+nuw1)/(V1_D+(Sigmas_D+nuw1));
    
    tau1_P=1/V1_P*(1-tildevz1_P/V1_P);
    tau1_D=1/V1_D*(1-tildevz1_D/V1_D);
    [tau1_P,tau1_P_old]=damping(tau1_P,tau1_P_old,mes);
    [tau1_D,tau1_D_old]=damping(tau1_D,tau1_D_old,mes);
    
    Sigmax_D=(N2*tau1_D*qh1)^(-1);
    Sigmah=(K*(rho*tau1_P+(1-rho)*qx1*tau1_D))^(-1);
    
    [Sigmax_D,Sigmax_D_old]=damping(Sigmax_D,Sigmax_D_old,mes1);
    [Sigmah,  Sigmah_old]=damping(Sigmah,Sigmah_old,mes1);
    
    
end

end