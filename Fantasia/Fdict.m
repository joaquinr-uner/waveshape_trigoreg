function [C] = Fdict(A,phi,D)
% Function that construct the pseudo-Fourier dictionary C given the 
% instantaneous amplitude (A) and phase (phi) and the number of
% harmonic components D.


C = zeros(length(A),2*D);

for i=1:D
    C(:,i) = A.*cos(2*pi*i*phi);
    C(:,D+i) = A.*sin(2*pi*i*phi);
end

end