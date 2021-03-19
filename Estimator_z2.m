function [tilde_z,tilde_vz]=Estimator_z2(Input,obj,Z,V,flag)
K1=Input.K(1);
if flag==1
    tilde_y=obj.tilde_Y(:,K1+1:end);
else   
    tilde_y=obj.tilde_Y;
end
nuw2=Input.nuw2;
[M,K]=size(tilde_y);
tilde_y=reshape(tilde_y,M*K,1);
V=reshape(V,M*K,1);
Z=reshape(Z,M*K,1);
if Input.ADC_switch==0
    tilde_vz=(V*nuw2)./(V+nuw2);
    tilde_z=(Z./V+tilde_y/nuw2)./(1/nuw2+1./V);
else %ADC case
    quan_step=obj.quan_step;
    Quan_bound=(2^(Input.n_bit-1)-1)*quan_step;
    tilde_y=[real(tilde_y);imag(tilde_y)];
    z=[real(Z); imag(Z)];
    v=[real(V); real(V)];
    y_up=tilde_y+quan_step/2;
    y_low=tilde_y-quan_step/2;
    [pos1,~]=find(y_up>Quan_bound);
    [pos2,~]=find(y_low<-Quan_bound);
    y_up(pos1)=1e5;
    y_low(pos2)=-1e5;
    
    eta1=(y_up-z)./sqrt((nuw2+v)/2);
    eta2=(y_low-z)./sqrt((nuw2+v)/2);
      
    tem1=normpdf(eta1)-normpdf(eta2);
    tem2=normcdf(eta1)-normcdf(eta2);
    tem3=eta1.*normpdf(eta1)-eta2.*normpdf(eta2);
    
    
    z_tem=z-v./sqrt(2*(nuw2+v)).*(tem1./tem2);
    v_tem=v/2-v.^2./(2*(nuw2+v)).*(tem3./tem2+(tem1./tem2).^2);
    
    tilde_z=z_tem(1:M*K)+1j*z_tem(M*K+1:end);
    tilde_vz=max(v_tem(1:M*K)+v_tem(M*K+1:end),eps);
    
end

tilde_z=reshape(tilde_z,M,K);
tilde_vz=reshape(tilde_vz,M,K); 

end

