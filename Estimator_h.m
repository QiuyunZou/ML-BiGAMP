function [hat_h,var_vh]=Estimator_h(Input,R_h,Sigma_h)

N1=Input.N1;
var_vh=Sigma_h./(1+N1*Sigma_h);
hat_h=(R_h./Sigma_h)./(1./Sigma_h+N1);

end