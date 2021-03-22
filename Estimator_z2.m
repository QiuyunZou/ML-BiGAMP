function [hatz, varz]=Estimator_z2(Input, obj, Z, V)
nuw2=Input.nuw2;
Y=obj.Y;
gain=V./(V+nuw2);
hatz=gain.*(Y-Z)+Z;
varz=nuw2*V./(V+nuw2);
end

