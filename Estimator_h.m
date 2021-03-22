function [hat_h,var_vh]=Estimator_h(R,Sigma)
var_vh=Sigma./(Sigma+1);
hat_h=R./(Sigma+1);
end