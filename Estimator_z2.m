function [tilde_z,tilde_vz]=Estimator_z2(Input,obj,Z,V)

tilde_y=obj.tilde_Y;
nuw2=Input.nuw2;
[N3,K]=size(tilde_y);
tilde_y=reshape(tilde_y,N3*K,1);
V=reshape(V,N3*K,1);
Z=reshape(Z,N3*K,1);


if Input.ADC_switch==0
    tilde_vz=(V*nuw2)./(V+nuw2);
    tilde_z=(Z./V+tilde_y/nuw2)./(1/nuw2+1./V);
else %ADC case
    quan_step=obj.quan_step;
    Quan_bound=(2^(Input.n_bit-1)-1)*quan_step;
    tilde_y=[real(tilde_y);imag(tilde_y)];
    z=[real(Z);imag(Z)];
    v=[0.5*real(V);0.5*real(V)];
    y_up=tilde_y+quan_step/2;
    y_low=tilde_y-quan_step/2;
    [pos1,~]=find(y_up>Quan_bound);
    [pos2,~]=find(y_low<-Quan_bound);
    y_up(pos1)=1e5;
    y_low(pos2)=-1e5;
    
    eta1=(sign(tilde_y).*z-min(abs(y_up),abs(y_low)))./sqrt(nuw2+v);
    eta2=(sign(tilde_y).*z-max(abs(y_up),abs(y_low)))./sqrt(nuw2+v);
      
    tem1=normpdf(eta1)-normpdf(eta2);
    tem2=normcdf(eta1)-normcdf(eta2);
    tem3=eta1.*normpdf(eta1)-eta2.*normpdf(eta2);
    
    pos=eta2<-100; 
    tem1(pos)=normpdf(eta1(pos));
    tem2(pos)=normcdf(eta1(pos));
    tem3(pos)=eta1(pos).*normpdf(eta1(pos));
    
    z_tem=z+(sign(tilde_y).*v./sqrt(nuw2+v)).*(tem1./tem2);
    v_tem=v-((v.^2)./(nuw2+v)).*(tem3./tem2+(tem1./tem2).^2);
    
    tilde_z=z_tem(1:N3*K)+1j*z_tem(N3*K+1:end);
    tilde_vz=max(v_tem(1:N3*K)+v_tem(N3*K+1:end),1e-10);
    
end

tilde_z=reshape(tilde_z,N3,K);
tilde_vz=reshape(tilde_vz,N3,K); 

end

