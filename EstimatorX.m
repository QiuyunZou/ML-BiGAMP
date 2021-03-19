function [hat_x,var_vx]=EstimatorX(obj, x_mean, x_var)
Allsymbol=obj.xo;
[N,K]=size(x_mean);
x_var=reshape(x_var,N*K,1);
x_mean=reshape(x_mean,N*K,1);

log_posterior=bsxfun(@times,-1./x_var,abs(bsxfun(@minus,x_mean,Allsymbol).^2));
log_posterior=bsxfun(@minus,log_posterior,max(log_posterior));  %防止溢出
posterior=exp(log_posterior); 
posterior=bsxfun(@rdivide,posterior,sum(posterior,2));           %得到标准PDF
hat_x=sum(bsxfun(@times,posterior,Allsymbol),2);                      %计算PDF的均值
var_vx=sum(posterior.*abs(bsxfun(@minus,hat_x,Allsymbol).^2),2);      %计算PDF的方差

hat_x=reshape(hat_x,N,K);
var_vx=reshape(var_vx,N,K);
end