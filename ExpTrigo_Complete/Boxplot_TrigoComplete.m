%% Generate the boxplot from Fig. 3 in Ruiz, J. Colominas, MA. 
% "Waveshape Model Order Estimation by Trigonometric Regression" 

res = 600;
SNR = {'0dB','5dB','10dB','15dB'};
vD0 = ['1','3','4','6','9'];
for i=1:length(vD0)
    load(['results\Exp_Trigo_Pure_' vD0(i) '.mat']);
    aux1(i) = aux;
    load(['results\Exp_Trigo_CLin_' vD0(i) '.mat']);
    aux2(i) = aux;
    load(['results\Exp_Trigo_CSin_' vD0(i) '.mat']);
    aux3(i) = aux;
end

figure(1);

ymin = min(min([aux1.DGCV,aux1.DRl,aux1.DWang,aux1.DKavalieris,aux2.DGCV,aux2.DRl,aux2.DWang,aux2.DKavalieris,aux3.DGCV,aux3.DRl,aux3.DWang,aux3.DKavalieris]));
ymax = max(max([aux1.DGCV,aux1.DRl,aux1.DWang,aux1.DKavalieris,aux2.DGCV,aux2.DRl,aux2.DWang,aux2.DKavalieris,aux3.DGCV,aux3.DRl,aux3.DWang,aux3.DKavalieris]));

DGCV = zeros(size(aux1(1).DGCV(1,:)));
DRl = DGCV;
DWang = DGCV;
DKavalieris = DGCV;

r_ds = -2.4;
primlabl = {'g', 'r', 'w','k'};
DBlabel = {'0dB','10dB','20dB','30dB'};
for j=1:size(aux1(1).DWang,1)
    for i=1:size(aux1,2)
        DGCV(i,:) = aux1(i).DGCV(j,:);
        DRl(i,:) = aux1(i).DRl(j,:);
        DWang(i,:) = aux1(i).DWang(j,:);
        DKavalieris(i,:) = aux1(i).DKavalieris(j,:);
    end
    h = {DGCV',DRl',DWang',DKavalieris'};
    
    subplot(4,3,3*(j-1)+1)
    h1 = boxplotGroup(h, 'PrimaryLabels', primlabl, ...
        'SecondaryLabels',{'1', '3', '4', '6', '9'}, 'InterGroupSpace', 1, 'OutlierSize', 3);
    set(gca,'YLim',[0 14])
    set(gca,'FontSize',12)
    ah = get(gca,'XLabel');
    set(ah, 'Position',[0 r_ds])
    xlabel('r_0','interpreter','tex')
    ah = get(gca,'YLabel');
    set(ah,'FontSize',24)
    if (j==1)
       title('No Modulation')
       ah = get(gca,'Title');
       set(ah,'FontSize',16)
    end
end

for j=1:size(aux1(1).DWang,1)
    for i=1:size(aux1,2)
        DGCV(i,:) = aux2(i).DGCV(j,:);
        DRl(i,:) = aux2(i).DRl(j,:);
        DWang(i,:) = aux2(i).DWang(j,:);
        DKavalieris(i,:) = aux2(i).DKavalieris(j,:);
    end
    h = {DGCV',DRl',DWang',DKavalieris'};
    subplot(4,3,3*(j-1)+2)
    h1 = boxplotGroup(h, 'PrimaryLabels', primlabl, ...
        'SecondaryLabels',{'1', '3', '4', '6', '9'}, 'InterGroupSpace', 1, 'OutlierSize', 3);
    set(gca,'YLim',[0 14])
    set(gca,'FontSize',12)
    ah = get(gca,'XLabel');
    set(ah, 'Position',[0 r_ds])
    xlabel('r_0','interpreter','tex')
    ah = get(gca,'YLabel');
    set(ah,'FontSize',24)
    
    if (j==1)
       title('Linear Modulation')
       ah = get(gca,'Title');
       set(ah,'FontSize',16)
    end
end

for j=1:size(aux1(1).DWang,1)
    for i=1:size(aux1,2)
        DGCV(i,:) = aux3(i).DGCV(j,:);
        DRl(i,:) = aux3(i).DRl(j,:);
        DWang(i,:) = aux3(i).DWang(j,:);
        DKavalieris(i,:) = aux3(i).DKavalieris(j,:);
    end
    h = {DGCV',DRl',DWang',DKavalieris'};
    subplot(4,3,3*(j-1)+3)
    h1 = boxplotGroup(h, 'PrimaryLabels', primlabl, ...
        'SecondaryLabels',{'1', '3', '4', '6', '9'}, 'InterGroupSpace', 1, 'OutlierSize', 3);
    set(gca,'YLim',[0 14])
    set(gca,'FontSize',12)
    ah = get(gca,'XLabel');
    set(ah, 'Position',[0 r_ds])
    xlabel('r_0','interpreter','tex')
    ah = get(gca,'YLabel');
    set(ah,'FontSize',24)
    if (j==1)
       title('Sinusoidal Modulation')
       ah = get(gca,'Title');
       set(ah,'FontSize',16)
    end
end