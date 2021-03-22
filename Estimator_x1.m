function [umean,uvar]=Estimator_x1(obj, m, v)
xo=obj.xo;
logpxr = bsxfun(@times, -1./v, abs(bsxfun(@minus, m, xo)).^2);
logpxr = bsxfun(@minus, logpxr, max(logpxr) );  % for stability          
pxr = exp(logpxr);
pxr = bsxfun(@rdivide, pxr, sum(pxr,2) );
umean = sum(bsxfun(@times, pxr, xo), 2); 
uvar = sum(pxr .* abs(bsxfun(@minus, umean, xo)).^2, 2);
end