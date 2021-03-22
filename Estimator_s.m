function [hat_sData, var_vsData]=Estimator_s(Input, obj, R_sData,Sigma_sData,Z1,V1)
nuw1=Input.nuw1;
omega=obj.omega;
var_vsData=1./(1./(nuw1+V1)+1./Sigma_sData);
hat_sData=var_vsData.*(Z1./(nuw1+V1)+R_sData./Sigma_sData);
var_vsData(~omega)=0;
hat_sData(~omega)=0;

end