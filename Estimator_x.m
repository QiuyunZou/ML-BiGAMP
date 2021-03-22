function [hatx, varx]=Estimator_x(R, Sigma)

hatx=R./(Sigma+1);
varx=Sigma./(Sigma+1);

end