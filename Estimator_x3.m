function [hat_x3, var_vx3]=Estimator_x3(Input,Z2, V2, R_x3, Sigma_x3)

nuw=Input.nuw;
var_vx3=1./(1./(nuw+V2)+1./Sigma_x3);
hat_x3=var_vx3.*(Z2./(nuw+V2)+R_x3./Sigma_x3);

end