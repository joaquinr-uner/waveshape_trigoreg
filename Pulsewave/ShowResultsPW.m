%% This file reproduces the boxplot from Fig. 7 from from the paper 
% Ruiz, J. Colominas, MA. "Waveshape Model Order Estimation by
% Trigonometric Regression".

    DB ={
     'AbdAorta'
};

T = ceil(size(DB,1)/4);

Suf = {'W','Rl','GCV','Sch'};

for j=1:size(DB,1)
    set(gcf, 'Position', [1 157 762 527]);
    load(['Results_WSF_Pulsewave_' DB{j} '.mat'])
    SNRout(:,3,:) = [];
    SNRout = circshift(SNRout,[0 2 0]);
    for i=1:size(SNRout,2)
        aux = SNRout(:,i,:);
        aux = reshape(aux,size(SNRout,1),size(SNRout,3));
        h{i} = aux;
    end
    subplot(3,4,1:8)
    boxplotGroup(h, 'PrimaryLabels',{'H','S','G','R','W','K','3','6','9','12'} , ...
        'SecondaryLabels',{'SNR_{in} = 0 dB' 'SNR_{in} = 10 dB' 'SNR_{in} = 20 dB'...
        'SNR_{in} = 30 dB'}, 'InterGroupSpace', 2)
    ylabel('SNR out (dB)')
    set(gca,'Position',[0.13 0.47 0.775 0.516])
    set(gca,'FontSize',fnts)
    

end