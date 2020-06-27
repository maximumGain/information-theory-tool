function H = entropy(px)
%This function computes the entropy of a random variable X with probability
%distribution pX(x).
%Entropy H(X) is the amount of uncertainty about a random variable X.
%H(X) = -sum_x pX(x)log2 pX(x)
%log base 2 entropy
px = px(:);
assert( all(px >= 0),'p(x) are not all non-negative');
assert( abs( 1 - sum(px) ) < 1E-10, 'px does not sum to 1');

%change 0 to 1E-10 avoids px*log2(px) = NaN
px(px < 1E-10) = 1E-10;

H = -sum( px .* (log2(px)));
end