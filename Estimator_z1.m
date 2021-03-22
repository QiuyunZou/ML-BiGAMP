function [tilde_z1, tilde_vz1]=Estimator_z1(Input,Z1,V1,R_x2,Sigma_x2)
nuw1=Input.nuw1;
tilde_vz1=(V1.*(Sigma_x2+nuw1))./(V1+(Sigma_x2+nuw1));
tilde_z1=tilde_vz1.*(R_x2./(Sigma_x2+nuw1)+Z1./V1);
end