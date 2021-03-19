function MSE =MLBiGAMP_SE(Input,obj)

N=Input.N;
M=Input.M;
IterNum=Input.IterNum ;
beta=M/N;
ADC_switch=Input.ADC_switch;
rho=Input.rho;

n_bit=Input.bit;
quan_step=obj.quan_step;
%% ADC setups
y_set=(-2^n_bit+1:2:2^n_bit-1)*quan_step/2;  %r的取值集合
y_up=y_set+quan_step/2;                      %r的上界
y_up(end)=Inf;
y_low=y_set-quan_step/2;                     %r的下界
y_low(1)=-Inf;

nuw=Input.nuw;

chi_x=1;
chi_h=1/M;
chi_z=N*chi_x*chi_h;
q_h=chi_h;
q_x=0;
V=N*chi_h*(chi_x-q_x)-1e-5;


Gaussian=@(x,m,v) 1./sqrt(2*pi*v).*exp(-(x-m).^2./(2*v));

Tem=1;
for ii=1:IterNum 
     if ADC_switch==1
        for Index=1:length(y_set)      
            alpha(Index)=integral(@(xi) (Gaussian(y_up(Index),sqrt(chi_z-V)*xi,nuw+V)-Gaussian(y_low(Index),sqrt(chi_z-V)*xi,nuw+V)).^2 ...
                      ./(normcdf((y_up(Index)-sqrt(chi_z-V)*xi)/(sqrt(nuw+V)))-normcdf((y_low(Index)-sqrt(chi_z-V)*xi)/(sqrt(nuw+V)))+eps).*Gaussian(xi,0,1),-Inf,Inf);
        end
        alpha=sum(alpha);
        mse_z=V-V^2*alpha;
    
        q_z=chi_z-mse_z;
        Sigma_x=N*q_h*(chi_x-q_x)^2/(beta*(q_z-N*q_x*q_h));
        
        
        MSE(ii)=1-rho/(1+rho*Sigma_x)*integral(@(z) z.^2./(rho+(1+rho)*sqrt((1+rho*Sigma_x)/(rho*Sigma_x))*exp(-z.^2./(2*rho*Sigma_x))+eps).*Gaussian(z,0,1),-Inf, Inf);
        
        q_x=chi_x-MSE(ii);
        V=N*chi_h*(chi_x-q_x);
     else 
         Sigma_x=nuw+1/beta*Tem;   %H=1/M；
         %Sigma_x=(nuw+Tem)/beta;   %H=1/N;
         MSE(ii)=1-rho/(1+rho*Sigma_x)*integral(@(z) z.^2./(rho+(1+rho)*sqrt((1+rho*Sigma_x)/(rho*Sigma_x))*exp(-z.^2./(2*rho*Sigma_x))+eps).*Gaussian(z,0,1),-Inf, Inf);
         Tem=MSE(ii);
     end
end

end