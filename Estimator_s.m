function [hat_sData, var_vsData]=Estimator_s(Input,R_sData,Sigma_sData,Z1,V1)
nuw1=Input.nuw1;
var_vsData=1./(1./(nuw1+V1)+1./Sigma_sData);
hat_sData=var_vsData.*(Z1./(nuw1+V1)+R_sData./Sigma_sData);
end