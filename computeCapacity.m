function computeCapacity(varargin)
%Function computeCapacity computes the capacity of a binary input channel,
%using pygx with "all" input distributions px.
%For a discrete memoryless channel pY|X(y|x), the “information” capacity C of a discrete memoryless channel is:
%C = max_pX(x) I(X;Y). The capacity is the maximization over all possible input distributions. 
%An optimal p∗X(x) is called the capacity-achieving input distribution: p∗X(x) = arg max_pX(x) I(X; Y)
%Capacity properties: 0<=C<=log|\mathcal{X}| and 0<=C<=log|\mathcal{Y}|. 
%log|\mathcal{X}|=log|\mathcal{Y}|=1 for binary case.
%space.
assert(nargin <= 1,'Too many input arguments');
if isempty(varargin)
    pygx=[...
        1 0 ; ...
        0 1];
else
    pygx = varargin;
end

%compute mutual information for all input distributions
p = linspace(eps, 1-eps,1001);%p = pX(1)
IXY = zeros(size(p));
for ii = 1:length(p)
    px =[ 1-p(ii) ; p(ii)];
    IXY(ii) = computeMutualInformation(pygx, px );
end

%Capacity is maximum of mutual information
[C,ind] = max(IXY);
pstar = p(ind);
fprintf('Capacity C = %g\n',C);
fprintf('Capacity−achieving input distribution p* = %g\n',pstar);
plot(p,IXY)
hold on
plot(pstar,C,'bo')
text(pstar+0.01,C-0.03,'Capacity')
axis([0 1 0 1])
xlabel('Input distribution p');
ylabel('Mutual Information I(X;Y)');
grid on