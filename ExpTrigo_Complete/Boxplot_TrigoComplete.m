%Date = '29_03_20';
Date = '28_05_20';
addpath('C:\Program Files\R2016a\toolbox\export_fig')
res = 600;
SNR = {'0dB','5dB','10dB','15dB'};
vD0 = ['1','3','4','6','9'];
for i=1:length(vD0)
    load(['Results\Exp_Trigo_Pure_' vD0(i) '.mat']);
    aux1(i) = aux;
    load(['Results\Exp_Trigo_CLin_' vD0(i) '.mat']);
    aux2(i) = aux;
    load(['Results\Exp_Trigo_CSin_' vD0(i) '.mat']);
    aux3(i) = aux;
end

figure(1);
%set(gcf, 'Position', get(0, 'Screensize'));
%set(gcf, 'Color','white')

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
    %ylabel(SNR{j})
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
%     set(gcf,'renderer','painters')
%     print(gcf,['C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Results_Pure' SNR{j} '-Square'],'-depsc',sprintf('-r%d',res))
%     saveas(gcf,['C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Results_Pure'  SNR{j} '-Square.fig'])
end

%saveas(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Results_Pure','epsc')
set(gcf,'renderer','painters')
%print(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Results_Pure','-depsc',sprintf('-r%d',res))
%saveas(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Results_Pure.fig')

%figure(2)
%set(gcf, 'Position', get(0, 'Screensize'));
%set(gcf, 'Color','white')
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
    %ylabel(SNR{j})
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
%     set(gcf,'renderer','painters')
%     print(gcf,['C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Results_LinMod' SNR{j} '-Square'],'-depsc',sprintf('-r%d',res))
%     saveas(gcf,['C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Results_LinMod'  SNR{j} '-Square.fig'])
end
%saveas(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Results_LinMod','epsc')
set(gcf,'renderer','painters')
%print(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Results_LinMod','-depsc',sprintf('-r%d',res))
%saveas(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Results_LinMod.fig')

%figure(3)
%set(gcf, 'Position', get(0, 'Screensize'));
%set(gcf, 'Color','white')
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
    %ylabel(SNR{j})
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
%      set(gcf,'renderer','painters')
%     print(gcf,['C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Results_SinMod' SNR{j} '-Square'],'-depsc',sprintf('-r%d',res))
%     saveas(gcf,['C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Results_SinMod'  SNR{j} '-Square.fig'])
end
%saveas(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Results_SinMod','epsc')
set(gcf,'renderer','painters')
%print(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Results_SinMod4x3','-depsc',sprintf('-r%d',res))
%saveas(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Results_SinMod4x3.fig')
