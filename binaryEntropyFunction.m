function h = binaryEntropyFunction(p)
%This function computes the entropy for binary random variable X
%of pX(x=0) = p and pX(x=1) = 1-p.
%It supports an input [p] of one number, or of a vector.
assert( all(p >= 0) && all(p <= 1),'p must be 0 <= p <= 1');
p(p<1E-10) = 1E-10;%To avoid NaN
p(1-p<1E-10)= 1- 1E-10;%To avoid NaN
h = - p .* log2(p) - (1-p) .* log2(1-p);
% h(p==0) = 0;
% h(p==1) = 0;
end