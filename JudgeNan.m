function judge=JudgeNan(A)
[M,N]=size(A);
A=reshape(A,M*N,1);
if sum(isnan(A))>0
    judge=1;
else 
    judge=0;
end
end
