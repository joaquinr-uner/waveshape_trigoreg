function h1 = spectro(F,t,f)
% Function that plots a spectrogram given the STFT of a signal
% and the time (t) and frequency (f) axis
if nargin<3
    imagesc(abs(F).^2), axis xy
else
    imagesc(t,f,abs(F).^2), axis xy
end
colormap(1-gray)

end

