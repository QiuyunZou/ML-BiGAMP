function MSE_error=MLBiGAMP(obj,Input)

IterNum=Input.IterNum;
N=Input.N;
M=Input.M;
H=obj.H;
mes=Input.mes;
MSE_error=zeros(IterNum,1);

%% array initialization
hat_v=ones(N,1);
hat_x=zeros(N,1);
hat_s=zeros(M,1);

sqrH=abs(H).^2;
sqrHt=sqrH';
Ht=H';
V_old=N/M*ones(M,1);
Z_old=obj.y;

MSE_old=1;

for ii=1:IterNum
    
    %% Output Nodes
    V=sqrH*hat_v;
    Z=H*hat_x-hat_s.*V;
    [Z,Z_old]=damping(Z,Z_old,mes);
    [V,V_old]=damping(V,V_old,mes);
    [hat_s,hat_tau]=Estimator_z(Z,V,obj,Input);
    
    %% Input Nodes
    Sigma=(sqrHt*hat_tau).^(-1);
    R=hat_x+Sigma.*(Ht*hat_s);  
    [hat_x,hat_v]=Estimator_x(Input, R, Sigma);
    MSE=norm(hat_x-obj.x)^2/(norm(obj.x)^2);
    
    
    MSE_error(ii,1)=MSE;
    MSE_old=MSE;
    
     
end


end




function [hatx,varx]=Estimator_x(Input,m,v)
rho=Input.rho;
sigma_X=1/rho;

Gau=@(x,a,v) 1./sqrt(2*pi*v).*exp(-1./(2*v).*abs(x-m).^2);

C=(rho*Gau(0,m,v+sigma_X))./...
    (rho*Gau(0,m,v+sigma_X)+(1-rho)*Gau(0,m,v));

hatx=C.*m./(1+rho*v);
Ex2=C.*(v./(1+rho*v)+(m./(1+rho*v)).^2);
varx=Ex2-abs(hatx).^2;
end


function [hat_s,hat_tau]=Estimator_z(Z,V,obj,Input)
nuw=Input.nuw;
ADC_switch=Input.ADC_switch;
y=obj.y;
M=Input.M;

if ADC_switch==0
    hat_s=(y-Z)./(V+nuw); 
    hat_tau=1./(V+nuw);
else
    quan_step=obj.quan_step;
    Quan_bound=(2^(Input.bit-1)-1)*quan_step;
    y_up=y+quan_step/2;
    y_low=y-quan_step/2;
    [pos1,~]=find(y_up>Quan_bound);
    [pos2,~]=find(y_low<-Quan_bound);
    y_up(pos1)=1e5;
    y_low(pos2)=-1e5;
    
    eta1=(y_up-Z)./sqrt(nuw+V);
    eta2=(y_low-Z)./sqrt(nuw+V);
    
      
    tem1=normpdf(eta1)-normpdf(eta2);
    tem2=normcdf(eta1)-normcdf(eta2);
    tem3=eta1.*normpdf(eta1)-eta2.*normpdf(eta2);
    
    pos=eta2<-100; 
    tem1(pos)=normpdf(eta1(pos));
    tem2(pos)=normcdf(eta1(pos));
    tem3(pos)=eta1(pos).*normpdf(eta1(pos));
    
    hatz=Z-V./sqrt(nuw+V).*(tem1./tem2);
    hatv=V-V.^2./(nuw+V).*(tem3./tem2+(tem1./tem2).^2);
    
    
    hat_s=(hatz-Z)./V;
    hat_tau=-(hatv-V)./(V.^2);
           
end
    


end