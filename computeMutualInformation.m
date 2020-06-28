function I = computeMutualInformation(varargin)
%method 1: Compute mutual information I(X; Y) using pY|X(y|x) and pX(x).
%          I = computeMutualInformation(pygx,px)
%method 2: Compute mutual information I(X; Y) using pX,Y(x,y).
%          I = computeMutualInformation(pxy)
%
%Let X and Y be jointly distributed random variables.
%The mutual information I(X;Y) between X and Y is: I(X; Y) = H(X) − H(X|Y)
%Property: I(X; Y) = H(X) + H(Y) − H(X, Y)
assert(nargin <= 2,'Too many input arguments.');
if isempty(varargin)
    pygx = [0.711 0.289; 0.289 0.711];
    px = [1/4 3/4];
elseif nargin == 2 %pygx and px
    assert( all(varargin{1}(:) >= 0),'Input are not all non-negative');
    assert( all(varargin{2}(:) >= 0),'Input are not all non-negative');
    pygx = varargin{1};
    assert( all(abs(1-sum(pygx,2))<1E-10), 'Row of pygx does not sum to 1');
    px = varargin{2};
    assert( abs( 1 - sum(px) ) < 1E-10, 'px does not sum to 1');
    
else %compute using pxy
    assert( all(varargin{1}(:) >= 0),'Input are not all non-negative');
    pxy = varargin{1};
    assert(abs( 1 - sum(pxy,'all') ) < 1E-10,'pxy does not sum to 1');
    py = sum(pxy,1);
    px = sum(pxy,2);
    HX = entropy(px);
    HY = entropy(py);
    HXY = entropy(pxy);
    I = HX + HY - HXY;
    return
    
end


%compute pxgy using inputs pygx and px
%first, find joint distribution pxy:
[X,Y] = size(pygx);
pxy = repmat(px(:),1,Y) .* pygx;
%second, find conditional distribution pxgy:
py = sum(pxy,1);
pxgy = pxy' ./ repmat(py,X,1)'; %transpose: rows sum to 1

%compute mutual information I = H(X) - H(X|Y)
HX = entropy(px);
HXgY = computeConditionalEntropy(pxgy,py);
I = HX - HXgY;
end