function [tilde_z2, tilde_vz2]=Estimator_z2(Input, Z2, V2, R_x3, Sigma_x3)
nuw=Input.nuw;
tilde_vz2=(V2.*(Sigma_x3+nuw))./(V2+(Sigma_x3+nuw));
tilde_z2=tilde_vz2.*(R_x3./(Sigma_x3+nuw)+Z2./V2);
end