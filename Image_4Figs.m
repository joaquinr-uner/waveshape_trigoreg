%% This file reproduces Fig. 4 from the paper 
% Ruiz, J. Colominas, MA. "Waveshape Model Order Estimation by
% Trigonometric Regression". 
Sufix = 'LinSin';
date = '12-Nov';
SNR = 0;
load(['Results_ExpTrigo_2Chirps_' date '-2020_RSR_' num2str(SNR) 'dB_' Sufix '.mat'])
DG = SA.DGCV;
DR = SA.DRl;
DW = SA.DWang;
DK = SA.DKavalieris;
r1max = max(DG(1,:));
r2max = max(DG(2,:));

[MG,LG,marg1G,marg2G] = MakeMat(DG);
[MR,LR,marg1R,marg2R] = MakeMat(DR);
[MW,LW,marg1W,marg2W] = MakeMat(DW);
[MK,LK,marg1K,marg2K] = MakeMat(DK);

cm = 1-gray;
txtfont = 14;
lbfnt = 22;
xlaps2 = -360;

figure(1)
%set(gcf,'Position',[0 171 560 513])
subplot(4,3,[1,2,4,5,7,8])
imagesc(LG(2,1):LG(2,2),LG(1,1):LG(1,2),MG), axis xy
xl = get(gca,'XLim');
yl = get(gca,'YLim');
posI = get(gca,'Position');
ylab=ylabel('$r^*_1$','interpreter','latex');
title(['GCV; ' num2str(SNR) ' dB'])
set(gca,'XTick',[]);set(gca,'YTick',LG(1,1):LG(1,2))
colormap(cm)
set(gca,'FontSize',txtfont)
subplot(4,3,[3,6,9])
bar(LG(1,1):LG(1,2),marg2G,'k')
posV = get(gca,'Position');
posV(2) = posI(2);
posV(4) = posI(4);
set(gca,'XLim',yl,'Position',posV)
camroll(-270)
set(gca,'FontSize',txtfont)
subplot(4,3,[10,11])
bar(LG(2,1):LG(2,2),marg1G,'k')
set(gca,'Position',[0.13 0.14 0.4942 0.1577])
set(gca,'XLim',xl)
xlab = xlabel('$r^*_2$','interpreter','latex');
xlab.Position(2) = xlaps2;
set(gca,'FontSize',txtfont)
set(ylab,'FontSize',lbfnt)
set(xlab,'FontSize',lbfnt)

figure(2)
%set(gcf,'Position',[0 171 560 513])
subplot(4,3,[1,2,4,5,7,8])
imagesc(LR(2,1):LR(2,2),LR(1,1):LR(1,2),MR), axis xy
xl = get(gca,'XLim');
yl = get(gca,'YLim');
posI = get(gca,'Position');
ylab=ylabel('$r^*_1$','interpreter','latex');
title(['Rl; ' num2str(SNR) ' dB'])
set(gca,'XTick',[]);set(gca,'YTick',LR(1,1):LR(1,2))
colormap(cm)
set(gca,'FontSize',txtfont)
subplot(4,3,[3,6,9])
bar(LR(1,1):LR(1,2),marg2R,'k')
posV = get(gca,'Position');
posV(2) = posI(2);
posV(4) = posI(4);
set(gca,'XLim',yl,'Position',posV)
camroll(-270)
set(gca,'FontSize',txtfont)
subplot(4,3,[10,11])
bar(LR(2,1):LR(2,2),marg1R,'k')
set(gca,'Position',[0.13 0.14 0.4942 0.1577])
set(gca,'XLim',xl)
xlab = xlabel('$r^*_2$','interpreter','latex');
xlab.Position(2) = xlaps2;
set(gca,'FontSize',txtfont)
set(ylab,'FontSize',lbfnt)
set(xlab,'FontSize',lbfnt)

figure(3)
%set(gcf,'Position',[0 171 560 513])
subplot(4,3,[1,2,4,5,7,8])
imagesc(LW(2,1):LW(2,2),LW(1,1):LW(1,2),MW), axis xy
xl = get(gca,'XLim');
yl = get(gca,'YLim');
posI = get(gca,'Position');
ylab=ylabel('$r^*_1$','interpreter','latex');
title(['Wang; ' num2str(SNR) ' dB'])
set(gca,'XTick',[]);set(gca,'YTick',LW(1,1):LW(1,2))
colormap(cm)
set(gca,'FontSize',txtfont)
subplot(4,3,[3,6,9])
bar(LW(1,1):LW(1,2),marg2W,'k')
posV = get(gca,'Position');
posV(2) = posI(2);
posV(4) = posI(4);
set(gca,'XLim',yl,'Position',posV)
camroll(-270)
set(gca,'FontSize',txtfont)
subplot(4,3,[10,11])
set(gca,'Position',[0.13 0.14 0.4942 0.1577])
bar(LW(2,1):LW(2,2),marg1W,'k')
set(gca,'XLim',xl)
xlab = xlabel('$r^*_2$','interpreter','latex');
xlab.Position(2) = xlaps2;
set(gca,'FontSize',txtfont)
set(ylab,'FontSize',lbfnt)
set(xlab,'FontSize',lbfnt)

figure(4)
subplot(4,3,[1,2,4,5,7,8])
imagesc(LK(2,1):LK(2,2),LK(1,1):LK(1,2),MK), axis xy
xl = get(gca,'XLim');
yl = get(gca,'YLim');
posI = get(gca,'Position');
ylab=ylabel('$r^*_1$','interpreter','latex');
title(['Kavalieris; ' num2str(SNR) ' dB'])
set(gca,'XTick',[]);set(gca,'YTick',LK(1,1):LK(1,2))
colormap(cm)
set(gca,'FontSize',txtfont)
subplot(4,3,[3,6,9])
bar(LK(1,1):LK(1,2),marg2K,'k')
posV = get(gca,'Position');
posV(2) = posI(2);
posV(4) = posI(4);
set(gca,'XLim',yl,'Position',posV)
camroll(-270)
set(gca,'FontSize',txtfont)
subplot(4,3,[10,11])
set(gca,'Position',[0.13 0.14 0.4942 0.1577])
bar(LK(2,1):LK(2,2),marg1K,'k')
set(gca,'XLim',xl)
xlab = xlabel('$r^*_2$','interpreter','latex');
xlab.Position(2) = xlaps2;
set(gca,'FontSize',txtfont)
set(ylab,'FontSize',lbfnt)
set(xlab,'FontSize',lbfnt)
