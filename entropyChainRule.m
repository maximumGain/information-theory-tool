function [HYgX,HXgY] = entropyChainRule(pxy)
%pxy, as a joint probability distribution pX,Y(x,y), is a matrix where the rows
%correspond to values of x and columns correspond to values of y.
%Chain Rule for Entropy For random variables X and Y:
%H(X, Y) = H(X) + H(Y|X) = H(Y) + H(X|Y)

assert( all(pxy(:) >= 0),'pxy are not all non-negative');
assert(abs( 1 - sum(pxy,'all') ) < 1E-10, 'pxy does not sum to 1');

%The marginal distribution of X and Y are obtained by summing columns and
%rows of pXY(x,y) respectively:
px = sum(pxy,2);
py = sum(pxy,1);
assert( abs( 1 - sum(px) ) < 1E-10, 'px does not sum to 1');
assert( abs( 1 - sum(py) ) < 1E-10, 'py does not sum to 1');

HX = entropy(px);
HY = entropy(py);
HXY = entropy(pxy);

HYgX = HXY - HX;
HXgY = HXY - HY;
end