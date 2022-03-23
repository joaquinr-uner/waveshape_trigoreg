addpath('C:\Program Files\R2016a\toolbox\export_fig')

%close all
    DB ={
     'AbdAorta';
};

labels = {'AA';'AT';'AR';'BR';'CA';'CI';'DI';'FE';'IB';'RA';'SC';'ST';'TA'};
c = 2.1;
T = ceil(size(DB,1)/4);

Suf = {'W','Rl','GCV','Sch'};
%Suf = 'Rl';

if(strcmp(Suf(1),'W')&&mod(c,1))
    cent = floor(c);
    cfrac = c - floor(c);
    cstg = [num2str(cent) 'dot' num2str(10*cfrac)];
else
    cstg = num2str(c);
end

%load(['C:\Users\Intel\Dropbox\Cluster\Exp_PW\Results_' DB{j} 'All-11-06-20r6.mat'])
%load(['C:\Users\Intel\Dropbox\Cluster\Exp_PW\Results_' DB{1} 'WithThr-17-07-20.mat'])
load(['E:\Resultados viejos Pulsewave Database\Results_' DB{1} 'WithThr-17-07-20.mat'])
names = {'H','S','G','R','W','K','3','6','9','12'};
SNRout = Meta.SNRout;
SNRout(:,3,:) = [];
SNRout = circshift(SNRout,[0 2 0]);
indcomb = combntns(1:10,2);
labcomb = combntns(names,2);

Mp = zeros(size(indcomb,1),size(SNRout,3));

alpha = 0.05;
alpha_bon = alpha/(2*size(indcomb,1));

SNR = [0,10,20,30];

MH = zeros(size(Mp));
Mdiff = zeros(size(Mp));
fileID = fopen('Test_Log.txt','w');
fprintf(fileID,'Hyphotesis H0: Pair Medians are equal. Hypothesis H1: Par Medians are unequal.\n');
fclose(fileID);
fileID = fopen('Test_Log.txt','a');
for j=1:size(SNRout,3)
    fprintf(['SNR: ' num2str(SNR(1)) 'dB.\n'])
    fprintf(fileID,['SNR: ' num2str(SNR(j)) 'dB.\n']);
for i=1:size(indcomb,1)
    Pair = SNRout(:,indcomb(i,:),j);
    %boxplot(Pair)
    Mdiff(i,j) = median(Pair(:,1))-median(Pair(:,2));
    Mp(i,j) = signrank(Pair(:,1),Pair(:,2));
    %Mp(i,j) = ranksum(Pair(:,1),Pair(:,2));
    
    if (Mp(i,j)<alpha_bon)
        MH(i,j) = 1;
    end
    fprintf([labcomb{i,1} '-' labcomb{i,2} '. p-value: ' num2str(Mp(i,j)) '. take Hypothesis ' num2str(MH(i,j)) ' Med.Diff:' num2str(Mdiff(i,j)) '\n'])
    fprintf(fileID,[labcomb{i,1} '-' labcomb{i,2} '. p-value: ' num2str(Mp(i,j)) '. Take Hypothesis H' num2str(MH(i,j)) ' Med.Diff:' num2str(Mdiff(i,j)) '\n']);
end
end

M0dB = cell(11,11);
M10dB = M0dB;
M20dB = M0dB;
M30dB = M0dB;

aux = {'$p-value$','H','S','$GCV$','$R_l$','Wang','Kavalieris','$r=3$','$r=6$','$r=9$','$r=12$'};

M0dB(1,:) = aux;
M10dB(1,:) = aux;
M20dB(1,:) = aux;
M30dB(1,:) = aux;

M0dB(:,1) = aux;
M10dB(:,1) = aux;
M20dB(:,1) = aux;
M30dB(:,1) = aux;
for i=1:size(indcomb,1)
   aux = Mp(i,1);
   if(aux<0.0001)
       
   M0dB{indcomb(i,2)+1,indcomb(i,1)+1} = '$<0.0001$';
   M0dB{indcomb(i,1)+1,indcomb(i,2)+1} = '$<0.0001$';
   else
   M0dB{indcomb(i,2)+1,indcomb(i,1)+1} = aux;
   M0dB{indcomb(i,1)+1,indcomb(i,2)+1} = aux;
   end
   
   aux = Mp(i,2);
   if(aux<0.0001)
       
   M10dB{indcomb(i,2)+1,indcomb(i,1)+1} = '$<0.0001$';
   M10dB{indcomb(i,1)+1,indcomb(i,2)+1} = '$<0.0001$';
   else
   M10dB{indcomb(i,2)+1,indcomb(i,1)+1} = aux;
   M10dB{indcomb(i,1)+1,indcomb(i,2)+1} = aux;
   end
   aux = Mp(i,3);
   if(aux<0.0001)
       
   M20dB{indcomb(i,2)+1,indcomb(i,1)+1} = '$<0.0001$';
   M20dB{indcomb(i,1)+1,indcomb(i,2)+1} = '$<0.0001$';
   else
   M20dB{indcomb(i,2)+1,indcomb(i,1)+1} = aux;
   M20dB{indcomb(i,1)+1,indcomb(i,2)+1} = aux;
   end
   aux = Mp(i,4);
   if(aux<0.0001)
       
   M30dB{indcomb(i,2)+1,indcomb(i,1)+1} = '$<0.0001$';
   M30dB{indcomb(i,1)+1,indcomb(i,2)+1} = '$<0.0001$';
   else
   M30dB{indcomb(i,2)+1,indcomb(i,1)+1} = aux;
   M30dB{indcomb(i,1)+1,indcomb(i,2)+1} = aux;
   end
end
txtsz = 8;
SNR0=SNRout(:,:,1);
SNR10=SNRout(:,:,2);
SNR20=SNRout(:,:,3);
SNR30=SNRout(:,:,4);
res = 600;

fig5=figure(5);
ax51 = subplot(3,4,9,'parent',fig5);
set(gca,'Position',[0.151, 0.1600, 0.1566, 0.2157])
set(ax51,'ytick',[])
set(ax51,'xtick',[])
ax52 = subplot(3,4,10,'parent',fig5);
set(gca,'Position',[0.3425, 0.1600, 0.1566, 0.2157])
set(ax52,'ytick',[])
set(ax52,'xtick',[])
ax53 = subplot(3,4,11,'parent',fig5);
set(gca,'Position',[0.535, 0.1600, 0.1566, 0.2157])
set(ax53,'ytick',[])
set(ax53,'xtick',[])
ax54 = subplot(3,4,12,'parent',fig5);
set(gca,'Position',[0.73, 0.1600, 0.1566, 0.2157])
set(ax54,'ytick',[])
set(ax54,'xtick',[])

fig1=figure(1);
set(gcf,'color','w'); 
[p,hbl,Stats] = friedman(SNR0,1,'off');
%c1 = multcompare(Stats,'CType','hsd');
[c1,~,comp1] = multcompare(Stats,'CType','bonferroni');
ax1=gca;
set(ax1,'ytick',[1:10],'yticklabel',fliplr(names))
ylabel(ax1,'Groups')
xlabel(ax1,'mean column rank')
title(ax1,'')
set(ax1,'FontSize',txtsz)
%aux = ylim;
%txt = text(0.1,aux(2)+2, 'TopLeft', 'HorizontalAlignment', 'left','VerticalAlignment', 'top');
%set(txt,'String','B','FontSize',30)
%set(gcf,'renderer','painters')
xlab = xlabel('mean column rank');
xlab.Position(2) = -1.7;
ylim([0 11])
axcp = copyobj(ax1,fig5);
set(axcp,'Position',get(ax51,'position'));
set(ax51,'YLim',[0,11])

%saveas(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Multcompare_0dB','epsc')
%p = get(gca,'Position');
%set(gca,'Position',[p(1) p(2)+0.1 p(3)-0.1 p(4)-0.2]);
%print(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Paper V1\Figures\Multcompare_0dB','-depsc',sprintf('-r%d',res))
%saveas(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Paper V1\Figures\Multcompare_0dB.fig')

figure(2)
[p,hbl,Stats] = friedman(SNR10,1,'off');
set(gcf,'color','w'); 
%c2 = multcompare(Stats,'CType','hsd');
c2 = multcompare(Stats,'CType','bonferroni');
ax2 = gca;
set(gca,'ytick',[])
%set(gca,'ytick',[1:10],'yticklabel',fliplr(names))
%ylabel('Groups')
xlabel('mean column rank')
title('')
set(gca,'FontSize',txtsz)
%aux = ylim;
%txt = text(0.1,aux(2), 'TopLeft', 'HorizontalAlignment', 'left','VerticalAlignment', 'top');
%set(txt,'String','C','FontSize',30)
xlab = xlabel('mean column rank');
xlab.Position(2) = -1.7;
ylim([0 11])
axcp = copyobj(ax2,fig5);
set(axcp,'Position',get(ax52,'position'));
set(ax52,'YLim',[0,11])

%saveas(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Multcompare_10dB','epsc')
%p = get(gca,'Position');
%set(gca,'Position',[p(1) p(2)+0.1 p(3)-0.1 p(4)-0.2]);
%print(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Paper V1\Figures\Multcompare_10dB','-depsc',sprintf('-r%d',res))
%saveas(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Paper V1\Figures\Multcompare_10dB.fig')

figure(3)
[p,hbl,Stats] = friedman(SNR20,1,'off');
set(gcf,'color','w'); 
%c3 = multcompare(Stats,'CType','hsd');
c3 = multcompare(Stats,'CType','bonferroni');
ax3 = gca;
set(gca,'yticklabel',[])
%set(gca,'ytick',[1:10],'yticklabel',fliplr(names))
%ylabel('Groups')
xlabel('mean column rank')
title('')
set(gca,'FontSize',txtsz)
aux = ylim;
%txt = text(0.1,aux(2), 'TopLeft', 'HorizontalAlignment', 'left','VerticalAlignment', 'top');
%set(txt,'String','D','FontSize',30)
set(gcf,'renderer','painters')
xlab = xlabel('mean column rank');
xlab.Position(2) = -1.7;
ylim([0 11])
axcp = copyobj(ax3,fig5);
set(axcp,'Position',get(ax53,'position'));
set(ax53,'YLim',[0,11])

%saveas(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Multcompare_20dB','epsc')
%p = get(gca,'Position');
%set(gca,'Position',[p(1) p(2)+0.1 p(3)-0.1 p(4)-0.2]);
%print(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Paper V1\Figures\Multcompare_20dB','-depsc',sprintf('-r%d',res))
%saveas(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Paper V1\Figures\Multcompare_20dB.fig')

figure(4)
[p,hbl,Stats] = friedman(SNR30,1,'off');
set(gcf,'color','w'); 
%c4 = multcompare(Stats,'CType','hsd');
c4 = multcompare(Stats,'CType','bonferroni');
ax4 = gca;
set(gca,'yticklabel',[])
%set(gca,'ytick',[1:10],'yticklabel',fliplr(names))
%ylabel('Groups')
xlabel('mean column rank')
title('')
set(gca,'FontSize',txtsz)
aux = ylim;
%txt = text(0.1,aux(2), 'TopLeft', 'HorizontalAlignment', 'left','VerticalAlignment', 'top');
%set(txt,'String','E','FontSize',30)
xlab = xlabel('mean column rank');
xlab.Position(2) = -1.7;
ylim([0 11])
set(gcf,'renderer','painters')
axcp = copyobj(ax4,fig5);
set(axcp,'Position',get(ax54,'position'));
set(ax54,'YLim',[0,11])

%saveas(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Figuras\Multcompare_30dB','epsc')
%p = get(gca,'Position');
%set(gca,'Position',[p(1) p(2)+0.1 p(3)-0.1 p(4)-0.2]); % Para que el texto aparezca cuando cambio el tamaño de la figure
%print(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Paper V1\Figures\Multcompare_30dB','-depsc',sprintf('-r%d',res))
%saveas(gcf,'C:\Users\Intel\Dropbox\ExpTrigo\Paper V1\Figures\Multcompare_30dB.fig')

%exl = actxserver('excel.application');
%set(exl.Selection,'HorizontalAlignment',3)
%xlswrite('M0dB.xlsx', M0dB);
%xlswrite('M10dB.xlsx', M10dB);  
%xlswrite('M20dB.xlsx', M20dB);  
%xlswrite('M30dB.xlsx', M30dB);  