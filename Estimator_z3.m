function [hatz, hatv]=Estimator_z3(Input, obj, m, v)
nuw=Input.nuw;
ADC_switch=Input.ADC_switch;
y=obj.tilde_Y;
N4=Input.N4;
DeltaTh=obj.DeltaTh;
quan_step=obj.quan_step;

if ADC_switch==0
    hatv=1./(1./v+1/nuw);
    hatz=hatv.*(y/nuw+m./v);
else
    y=[real(y);imag(y)];
    m=[real(m);imag(m)];
    v=[real(v);real(v)];
    y_up=y+quan_step/2;
    y_low=y-quan_step/2;
    [pos1,~]=find(y>max(DeltaTh));   
    [pos2,~]=find(y<-max(DeltaTh));
    y_up(pos1)=1e5;
    y_low(pos2)=-1e5;

    
    eta1=(y_up-m)./sqrt((nuw+v)/2);
    eta2=(y_low-m)./sqrt((nuw+v)/2);
      
    tem1=normpdf(eta1)-normpdf(eta2);
    tem2=normcdf(eta1)-normcdf(eta2);
    tem3=eta1.*normpdf(eta1)-eta2.*normpdf(eta2);
    
    
    z_tem=m-v./(sqrt(2*(nuw+v))).*(tem1./tem2);
    v_tem=v/2-((v.^2)./(2*(nuw+v))).*(tem3./tem2+(tem1./tem2).^2);
    
    hatz=z_tem(1:N4)+1j*z_tem(N4+1:2*N4);
    hatv=max(v_tem(1:N4)+v_tem(N4+1:2*N4),eps);         
end
end