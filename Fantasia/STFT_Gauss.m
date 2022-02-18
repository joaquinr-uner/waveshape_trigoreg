function [ F, sG ] = STFT_Gauss( x,K,sigma,fmax )
%STFT_Gauss implements the Short-Time Fourier Transform (STFT) using a
%Gaussian window in the frequency domain
X = fft(x);
N = length(x);
k = 0:1/K:fmax-1/K;
f = 0:1/N:1-1/N;
sG = 0;
for i=1:length(k)
    G = sqrt(pi/sigma)*exp(-pi^2*(f-k(i)).^2/sigma);
    G = G';
    F(i,:) = ifft(X.*G);
    sG = sG + abs(G);
end
end