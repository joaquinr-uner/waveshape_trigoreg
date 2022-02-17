%% This code reproduces Fig. 9 from the paper 
%% Ruiz, J. Colominas, MA. "Waveshape Model Order Estimation by
%% Trigonometric Regression" and the Figures from the 
%% Supplemental Material of said paper.

drt_r = 'results';
drt_d = 'data';

listing_r = dir(drt_r);
listing_d = dir(drt_d);
vst=[repmat(1e5,[1,8]),2.5e5,repmat(1e5,[1,11]),repmat(7.5e5,[1,2]),1e5,1,...
    1e5,7.1e5,repmat(1e5,[1,2]),7.7e5,repmat(1e5,[1,4]),7.1e5,repmat(1e5,[1,3]),...
    7.5e5,1e5,7.5e5];
ved=[repmat(1.5e5-1,[1,8]),3e5-1,repmat(1.5e5-1,[1,11]),repmat(8e5-1,[1,2]),...
    1.5e5-1,5e4,1.5e5-1,7.6e5-1,repmat(1.5e5-1,[1,2]),8.2e5-1,...
    repmat(1.5e5-1,[1,4]),7.6e5-1,repmat(1.5e5-1,[1,3]),8e5-1,1.5e5-1,8e5-1];

fs=250;
N=1e7;
t=0:1/fs:N/fs-1/fs;
res=600;
pos=[425 301 543 365];
for i=3:length(listing_d)
   dname = listing_d(i).name;
   rname = ['Results_WSF_' listing_d(i).name];
   load([drt_r '/' rname]);
   data = load([drt_d '/' dname]);
   st=vst(i-2);
   ed=ved(i-2);
   ecg = data.val(2,st:ed);
   resp = data.val(1,st:ed);
   ta=t(st:ed);
   tpi=linspace(-pi,pi,length(Rr.ECG_WSF));
   figure(1)
    set(gcf,'render','painters')
    set(gcf,'Position',pos);
   subplot(223)
   plot(tpi,Rr.ECG_WSF,'k')
   xlim([-pi,pi])
   ylim([min(Rr.ECG_WSF),max(Rr.ECG_WSF)])
   title(['Est. Order: ' num2str(length(Rr.ECG_Coeff)/2) '. Max. Order: ' num2str(Rr.rmax)])
   subplot(224)
   plot(tpi,Rr.Resp_WSF,'k')
   xlim([-pi,pi])
   ylim([min(Rr.Resp_WSF),max(Rr.Resp_WSF)])
   title(['Est. Order: ' num2str(length(Rr.Resp_Coef)/2) '. Max. Order: ' num2str(Rr.rmax_resp)])
   h3=subplot(221);
   %plot(ecg(1:1000)-mean(ecg(1:1000)))
   plot(ta,ecg-mean(ecg))
   hold on
   plot(ta,Rr.ECG_Estimate,'k')
   hold off
   set(gca,'YTick',[])
   xlim([ta(25002) ta(26002)])
   xlabel('Time (s)')
   ylabel('Amplitude')
   %legend('Original','Estimated')
   title(dname(1:end-4))
   h4=subplot(222);
   %plot(resp(1:1000)-mean(resp(1:1000)))
   plot(ta,resp-mean(resp))
   hold on
   plot(ta,Rr.Resp_Estimate,'k')
   hold off
   set(gca,'YTick',[])
   xlim([ta(25002) ta(32502)])
   xlabel('Time (s)')
   ylabel('Amplitude') 
   %legend('Original','Estimated')
   title(dname(1:end-4))
   %linkaxes([h3 h4],'x')
   print(gcf,['figures\Fig_' dname(1:end-4)],'-depsc',sprintf('-r%d',res))
end
