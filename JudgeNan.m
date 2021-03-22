function judge=JudgeNan(A)
[M,N]=size(A);
Ra=real(A);
Ia=imag(A);
Ra=reshape(Ra,M*N,1);
Ia=reshape(Ia,M*N,1);

Flag_R=sum(isnan(Ra)+isinf(Ra));
Flag_I=sum(isnan(Ia)+isinf(Ia));
Flag=Flag_R+Flag_I;
if Flag>0 
    judge=1;
else 
    judge=0;
end
end