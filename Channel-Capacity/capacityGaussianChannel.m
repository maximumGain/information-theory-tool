%This source code computes the capacity of AWGN channel, and shows the
%capacity with  binary input.
snrdB = -15:15; %snr in dB
VAR = 1./(10.^(snrdB./10));

biAWGN_C = zeros(length(VAR),1);
for ii = 1:length(VAR)
    var = VAR(ii);
    fun = @(y) (1/sqrt(2*pi*var)) .* exp(-(y-1).^2 / (2*var) ).* log2( 2 ./ ( 1 + exp(-2*y / var)) );
    biAWGN_C(ii) = integral(fun,-10,10);
end

snr = 1./VAR;
%Gaussian Channel Theorem P=1
AWGN_C = (1/2)*log2(1+snr);

figure
plot(snrdB,biAWGN_C,'b')
hold on
plot(snrdB,AWGN_C,'r')
grid on
xlabel('SNR in dB')
ylabel('Capacity')
legend(' binary-input AWGN channel','AWGN channel')
title('Capacity of AWGN Channel')