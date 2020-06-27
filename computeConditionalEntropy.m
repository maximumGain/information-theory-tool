function [HYgX,HYgXeqx] = computeConditionalEntropy(pygx,px)
%The conditional entropy H(Y|X) of discrete random variables X and Y jointly distributed as pX,Y(x,y) is:
% H(Y|X) = sum_x∈X  pX(x)H(Y|X = x)
%pygx is a matrix of length(px) rows and length(py) columns.
%The conditional entropy H(Y|X = x) is given by:
% H(Y|X = x) = -sum_y∈Y pY|X(y|x) log pY|X(y|x)
assert( all(px >= 0),'px are not all non-negative');
assert( abs( 1 - sum(px) ) < 1E-10, 'px does not sum to 1');
%The conditional entropy is usually written in a form so that each row sum
%to 1.
assert( all(abs(1-sum(pygx,2))<1E-10), 'Row of pygx does not sum to 1');

[X,~] = size(pygx);
assert(X == length(px),'number X elements in px and pygx disagree')

HYgXeqx = zeros(1,X); %HYgXeqx(x) is H(Y | X = x)
for x = 1:X
HYgXeqx(x) = entropy(pygx(x,:));
end

HYgX = HYgXeqx * px(:); %H(Y|X) = sum_x p(x)*HYgXeqx(x)
end