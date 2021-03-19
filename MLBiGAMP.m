function MSE_error=MLBiGAMP(obj,Input)

IterNum=Input.IterNum;
N=Input.N;
M=Input.M;
H=obj.H;
xo=obj.xo;
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
Z_old=zeros(M,1);

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
    [hat_x,hat_v]=Estimator_x(xo,R, Sigma);
    MSE=norm(hat_x-obj.x)^2/N;
    Copy_MSE=mean(hat_v);
    
    if sum(isnan(hat_v))>0 || MSE>MSE_old
       MSE_error(ii:IterNum,1)=MSE_old;
       break;
    end
    
    MSE_error(ii,1)=MSE;
    MSE_old=MSE;
    
     
end


end




function [hatx,varx]=Estimator_x(xo, m, v)
log_posterior=bsxfun(@times,-1./v,abs(bsxfun(@minus,xo,m).^2));
log_posterior=bsxfun(@minus,log_posterior,max(log_posterior));  %防止溢出

posterior=exp(log_posterior); 
posterior=bsxfun(@rdivide,posterior,sum(posterior,2));       %得到标准PDF
hatx=sum(bsxfun(@times,posterior,xo),2);                     %计算PDF的均值
varx=sum(posterior.*abs(bsxfun(@minus,hatx,xo).^2),2);       %计算PDF的方差
end


function [hat_s,hat_tau]=Estimator_z(Z,V,obj,Input)
nuw=Input.nuw;
AGC_switch=Input.AGC_switch;
y=obj.y;
M=Input.M;

if AGC_switch==0
    hat_s=(y-Z)./(V+nuw); 
    hat_tau=1./(V+nuw);
else
    quan_step=obj.quan_step;
    Quan_bound=(2^(Input.bit-1)-1)*quan_step;
    y=[real(y);imag(y)];
    z=[real(Z);imag(Z)];
    v=[real(V); real(V)];
    y_up=y+quan_step/2;
    y_low=y-quan_step/2;
    [pos1,~]=find(y_up>Quan_bound);
    [pos2,~]=find(y_low<-Quan_bound);
    y_up(pos1)=1e5;
    y_low(pos2)=-1e5;
    
    eta1=(y_up-z)./sqrt((nuw+v)./2);
    eta2=(y_low-z)./sqrt((nuw+v)./2);
    
      
    tem1=normpdf(eta1)-normpdf(eta2);
    tem2=normcdf(eta1)-normcdf(eta2);
    tem3=eta1.*normpdf(eta1)-eta2.*normpdf(eta2);
    
    pos=eta2<-100; 
    tem1(pos)=normpdf(eta1(pos));
    tem2(pos)=normcdf(eta1(pos));
    tem3(pos)=eta1(pos).*normpdf(eta1(pos));
    
    z_tem=z-v./sqrt(2*(nuw+v)).*(tem1./tem2);
    v_tem=v/2-v.^2./(2*(nuw+v)).*(tem3./tem2+(tem1./tem2).^2);
    
    hatz=z_tem(1:M)+1j*z_tem(M+1:2*M);
    hatv=max(v_tem(1:M)+v_tem(M+1:2*M),1e-10);
    
    hat_s=(hatz-Z)./V;
    hat_tau=-(hatv-V)./(V.^2);
           
end
    


end