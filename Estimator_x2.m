function [hat_x2, var_vx2]=Estimator_x2(Input,Z1, V1, R_x2, Sigma_x2)

nuw=Input.nuw;
var_vx2=1./(1./(nuw+V1)+1./Sigma_x2);
hat_x2=var_vx2.*(Z1./(nuw+V1)+R_x2./Sigma_x2);
end