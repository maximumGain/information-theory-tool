function [RD,D] = rateDistortionComputation(px,d,lambda)
%Arimoto-Blahut algorithm in test channel
%x as input, xhat(y) as output
%px:probability distribution of x
%d:distortion measure d(x,xhat)
%D:the expected distortion for a (2^nR,n) code 
%initial random input distribution
%y=xhat=r
[X,Y] = size(d);
r = rand(1,Y);
r = r / sum(r);


R = 2*r; %initiate any R greater than r.
SYM = 1; %a symbol to control the loop.
q = zeros(X,Y);
while (SYM)
     
%Fix r(xhat), minimize over q(xhat|x):
for y = 1:Y     
for x = 1:X
    %q(y|x)
q(x,y) = r(y) * exp(-lambda*d(x,y)) / sum(r .*exp(-lambda*d(x,:)));
end
end
%Fix q(xhat|x), minimize over r(xhat):
r = sum(repmat(px',1,Y).*q);

if abs(R-r)>= 0.000001
   SYM = 1;
   R = r;
else
   SYM = 0;
end
end

qstar = q;
rstar = r;
%compute RateDistortion
RD = sum(sum(repmat(px',1,Y).*qstar.*log2(qstar./repmat(rstar,X,1))));
D = sum(sum(repmat(px',1,y).*qstar.*d));
return