function [tilde_z1, tilde_vz1]=Estimator_z1(Input,Z1,V1,R_s,Sigma_s)
nuw=Input.nuw;
tilde_vz1=(V1.*(Sigma_s+nuw))./(V1+(Sigma_s+nuw));
tilde_z1=tilde_vz1.*(R_s./(Sigma_s+nuw)+Z1./V1);
end