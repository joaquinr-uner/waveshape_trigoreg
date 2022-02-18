function [x_r, v] = WSF(x,c,D,mode,b,F,sF)

global sigma
global fmax
if nargin<6
    switch mode
        case 1
            [F,sF] = STFT_Gauss(x,length(x),sigma,fmax);
        case 2
            [F1, sF] = STFT_Gauss(x,length(x),sigma,fmax);
            [F] = synchro_F(F1,0.5);
        case 3
            [F1, sF] = STFT_Gauss(x,length(x),sigma,fmax);
            F = threshold(F1);
        otherwise
    end
end
y = zeros(size(x));

for i=1:length(y)
   binf = max([1,c(i)-b]);
   bup = min([size(F,1),c(i)+b]);
   y(i) = 1/max(sF)*sum(F(binf:bup,i)); 
end

A = abs(y);
phi=unwrap(angle(y))/(2*pi);


C = Fdict(A,phi,D);

v = ((C'*C)\C')*x;

x_r = C*v;
end