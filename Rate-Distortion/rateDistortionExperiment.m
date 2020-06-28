%Rate distortion experiment
px = [1/3 1/3 1/3];
d = [0 1 4;1 0 1;4 1 0];

lambda = linspace(10^-2,15,100);

RD = zeros(length(lambda),1);
D = zeros(length(lambda),1);
for ii = 1:length(lambda)
[RD(ii),D(ii)] = rateDistortionComputation(px,d,lambda(ii));
end

figure
plot(D,RD,'r')
axis([0 1 0 2])
grid on
xlabel('D')
ylabel('RD')
title('Rate Distortion vs D')

