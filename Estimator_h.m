function [hat_h,var_vh]=Estimator_h(Input,R_h,Sigma_h)
N=Input.N;
var_vh=1./(1./Sigma_h+N);
hat_h=(R_h./Sigma_h)./(1./Sigma_h+N);
end