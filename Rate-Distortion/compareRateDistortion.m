%rate distortion for Bernoulli p source
pmatrix = [0.1 0.2 0.5];
RD = zeros(length(pmatrix),length(lambda));
D = zeros(length(pmatrix),length(lambda));
% hD = zeros(length(pmatrix),length(lambda));
% RDhat = zeros(length(pmatrix),length(lambda));
for jj = 1:length(pmatrix)
    p = pmatrix(jj);
    d = [0 1;1 0];%The binary Hamming distortion function
    lambda = linspace(10^-2,15,100);
    %     hp = binaryEntropyFunction(p);
    
    for ii = 1:length(lambda)
        %     1: Arimoto-Blahut algorithm in test channel
        [RD(jj,ii),D(jj,ii)] = rateDistortionComputation([1-p p],d,lambda(ii));
        %     2: RDhat: Rate-distortion function for binary source
        %     hD(jj,ii) = binaryEntropyFunction(D(jj,ii));
        %     RDhat(jj,ii) = hp - hD(jj,ii);
    end
end
%RDhat from the formula
%hold on
%plot(D,RDhat,'bo')

%Rate distortion is a function of the expected distortion D
figure
plot(D(1,:), RD(1,:),'r-')
hold on
hold on
plot(D(2,:), RD(2,:),'g-')
hold on
plot(D(3,:), RD(3,:),'b-')
axis([0 0.6 0 1])
grid on
xlabel('D')
ylabel('RD')
legend('p = 0.1','p = 0.2','p = 0.5')
title('Rate−Distortion for Bernoulli p source')


