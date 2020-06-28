function [C,pxstar] = capacityComputationConvergence(varargin)
%This function shows convergence of computing capacity.
%It implements the Arimoto-Blahut algorithm for computing the capacity 
%of a DMC given by pY|X(y|x)
assert(nargin <= 1,'Too many input arguments');
if isempty(varargin)
    pygx=[...
        1 0 ; ...
        0 1];
else
    pygx = varargin;
end

assert( all(pygx(:) >= 0),'pygx are not all non-negative');
assert( all(abs(1-sum(pygx,2))<1E-10), 'Row of pygx does not sum to 1');

[X,Y] = size(pygx);

%initial random input distribution
px = rand(1,X);
px = px / sum(px);

pygx(pygx==0) = 1E-10;

r = px;
p = pygx;

R = 2*r; %initiate any R greater than r.
SYM = 1; %a symbol to control the loop.
q = zeros(X,Y);
while (SYM)
    
    %Fix r(x), maximize over q(x|y):
    for x = 1:X
        for y = 1:Y
            q(x,y) = p(x,y) * r(x) / sum(r(:) .* p(:,y));
        end
    end
    
    %Fix q(x|y), maximize over r(x):
    for x = 1:X
        r(x) = prod(q(x,:).^p(x,:));
    end
    r = r ./ sum(r);
    
    %check the convergence
    if abs(R-r)>= 0.000001
        SYM = 1;
        R = r;
    else
        SYM = 0;
    end
end

pxstar = r;
pxy = repmat(pxstar',1,Y) .* pygx; %joint distribution pxy
py = sum(pxy,1); %marginal distribution py

%Capacity is mutual information using pxstar
C = sum(sum(pxy .* log2( pxy ./ (pxstar' * py))));
end