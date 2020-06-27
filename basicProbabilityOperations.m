%--------------------------------1
%Jointly distributed random variables X and Y are defined for sample spaces X and Y respectively, 
%with a joint probability distribution: pX,Y(x, y) = Pr(X = x, Y = y)
%The joint distribution pX,Y (x, y) is the “all-knowing distribution.” 
%This means that given pX,Y(x,y), it is possible to compute pX(x), pY(y), pX|Y(x|y) and pY|X(y|x). 
%If you are given only pY|X(y|x), then you additionally need pX(x) to find
%the joint distribution.
%--------------------------------1

%"Theorem of total probability" PX(x) = sum_y pY(y)PX|Y(x|y) 
%compute PY(y) = sum_x pX(x)PY|X(y|x) 
%pygx is "probability of Y given X", as a matrix of size(X) rows and
%size(Y) columns. Rows sum to 1.
pygx = [...
          0     1     0;
        1/2     0   1/2;
        1/3   1/3   1/3;
        0     2/3   1/3];
assert(all(abs(1-sum(pygx,2))<1E-10),'Each row of pygx should sum to 1');
px = [1/4 3/8 1/8 1/4];
assert(abs(1-sum(px))<1E-10,'px does not sum to 1');
py = px * pygx;

%"Marginalization" PY(y) = sum_y PX,Y(x,y) and PX(x) = sum_x PX,Y(x,y)
pxy = [1/2 0 0; 1/4 1/16 0; 0 1/8 1/16];
py = sum(pxy,1);
px = sum(pxy,2);
rats(px) %display as rational numbers



%"Compute joint distribution" PX,Y(x,y) = PY|X(y,x)PX(x)
[X,Y] = size(pygx);
pxy = repmat(px(:),1,Y) .* pygx;

%"Compute the conditional probability" PX|Y(x|y) = PX,Y(x,y)/PY(y)
[X,Y] = size(pxy);
py = sum(pxy,1);
px = sum(pxy,2);
pxgy = pxy' ./ repmat(py(:)',X,1)';%transpose: rows sum to 1
%pxgy is a matrix of Y rows and X columns.
pygx = pxy ./ repmat(px(:),1,Y);

%"Bayes rule" PX|Y(x|y) = PY|X(y|x)PX(x)/PY(y)
pxgy = transpose(pygx .* repmat(px(:),1,Y) ./ repmat(py(:)',X,1));

%--------------------------------2
%Two random variables X and Y are independent if and only if: pX,Y(x,y) = pX(x)pY(y)
%Also, X and Y are independent if an only if pX|Y(x|y) = pX(x) for all x ∈ X,y ∈ Y.
%If pX(x) = pY(x) for all x ∈ X, then we say X and Y are independent and identically distributed, often abbreviated iid.
%--------------------------------2


%--------------------------------3
%Expectation and Variance
%--------------------------------3

%"Compute expected value" E[X] = sum_x x pX(x)
%die roll random variable X
x = 1:6; %X = {1,2,3,4,5,6}
px = repmat(1/6,6,1); %pX(x) = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6] 
EX = X*px(:);
%"Compute expected value of a function" g(x)] E[g(X)] = sum_x g(x)pX(x)
g = @(x) x.^2; %g(x) = x^2
gp = @(px) -log2(px); %gp(x) = -log pX(x)
x = [0 1 2 3];
px = [1/4 3/8 1/8 1/4];
EgX = g(x) * px(:);

%"compute variance" Var[X] = E[X^2] - (E[X])^2
%die roll random variable X
x = 1:6; %X = {1,2,3,4,5,6}
px = repmat(1/6,6,1); %pX(x) = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6] 
EX = x*px(:);
g = @(x) x.^2; %g(x) = x^2
EgX = g(x) * px(:);
Var = EgX - (EX)^2;

%"The conditional expectation" E[X|Y = y] = sum_x xPX|Y(x|y)
%The conditional expectation E[X|Y = y] is a function of y, for example, f(y) = E[X|Y = y]. 
%Note that E[X|Y] is distinct from E[X|Y = y].
%In particular E[X|Y] is a random variable equal to f(Y). 
%Since E[X|Y] is a random variable, it has an expectation. In fact, it is equal to E[X].
%Proposition: [Law of Total Expectation] Let X and Y be jointly distributed random variables
%E[E[X|Y]] = E[X]


%"Expectation of Sums of Random Variables"
%For any X1,X2,...,Xn and constants a1, a2,...,an, E[a1X1 + a2X2 + · · · anXn] = a1E[X1] + a2E[X2] + · · · anE[Xn]

%"Variance of Sums of Random Variables"
%For any independent X1,X2,...,Xn and constants a1, a2,...,an, Var[a1X1 + a2X2 + · · · anXn] = a21Var[X1] + a2Var[X2] + · · · a2nVar[Xn]



%--------------------------------4
%Random Vectors
%--------------------------------4
%Let X = X1,X2,...,Xn be a random vector of n random variables, independent and identically distributed.
%A random variable Xi has distribution pX(x). 
%Then, the random vector X =(X1,X2,X3,...,Xn)

%"The sample mean Xnbar" Xnbar = (1/n) sum Xi for i = 1,2,...,n
%"The expected mean of the binary random vector" E[Xn_bar] = E[Xi]

%"Binary random vector" If Xi is a binary random variable on sample space {0, 1}, 
%with probability of a one equal to 0 ≤ p ≤ 1, then we say X is a binary random vector.
%The expected mean of binary random vector X is E[Xn_bar] = E[Xi] = p.
%The variance of binary random vector X is Var[Xn_bar] = p(1-p)/n.

%The binary random vector X has k ones and n − k zeros. 
%Let K be the sum of X: K = sum Xi for i = 1,2,...,n, so that K is a random
%variable expressing the number of ones in X. This is the binomial random 
%variable K with probability distribution 
%pK(k) = (n!/k!(n-k)!)p^k (1-p)^(n-k) for k = 0,1,...,n


%"The General Binomial Probability Formula" P(k out of n) = pK(k) = (n!/k!(n-k)!) p^k (1-p)^(n-k)
%Great explanation at https://www.mathsisfun.com/data/binomial-distribution.html
pk = nchoosek(n,k) * p^k * (1-p)^(n-k);

%Example: consider an example of a binary random vector with n = 4 and p = 1/4 (pX(1) = p = 1/4 and pX(0) = 1-p = 3/4). 
%The sample space is: {0000,0001,0010,...,1111}.
%The probability of X = {0100} is pX0100 = pX(0)pX(1)pX(0)pX(1) by dependency.
p = 1/4;
pX0100 = (1-p)*p*(1-p)*(1-p);
%The probability that X is a sequence of two 1's and two 0's: 27/128
n = 4;
k = 2;
pk = nchoosek(n,k) * p^k * (1-p)^(n-k);
%The probability that X has one 1 or X has two 1's: 81/128
k = 1;
pk1 = nchoosek(n,k) * p^k * (1-p)^(n-k);
k = 2;
pk2 = nchoosek(n,k) * p^k * (1-p)^(n-k);
pk1or2 = pk1+pk2;








