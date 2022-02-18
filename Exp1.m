%% This file reproduces the experiment presented in Fig. 1 of 
% Ruiz, J. Colominas, MA. "Waveshape Model Order Estimation by 
% Trigonometric Regression". 
amp = [1,0.7,0.6,0.5,0.4,0.3,0.4,0.3,0.2];

res = 600;
vD0 = [4,9];

vRSR = [Inf,10,0];

N= 3000;
rng(1)
t = linspace(0,1,N);
vc = [1,2,5,8,12];
H = floor((log(N))^2);
fs = 1/(t(2)-t(1));

f = 0:fs/N:fs/2-fs/N;
G = sqrt(sqrt(pi)/sqrt(2*1e-4));

phi = 70*t + 15/(2*pi)*cos(2*pi*t);

Trig = 'Sin';

s4 = Fdict(ones(1,N),linspace(-1/2,1/2,N),4)*[amp(1:4) zeros(1,4)]';
s9 = Fdict(ones(1,N),linspace(-1/2,1/2,N),9)*[amp(1:9) zeros(1,9)]';

S4 = abs(fftshift(fft(s4)));
S4 = S4/max(S4);
S9 = abs(fftshift(fft(s9)));
S9 = S9/max(S9);

txtfont = 7;
figure(1)
set(gcf,'Position',[10 338 418 346])
subplot(221);plot(linspace(-pi,pi,length(s4)),s4,'k');xlim([-pi,pi])
xlabel('$t$','interpreter','latex')
ylim([-1 3])
set(gca,'FontSize',txtfont,'FontUnits','points')
text(-3,2.5,'$s_1(t)$','interpreter','latex','FontSize',12)

subplot(222);plot(linspace(-pi,pi,length(s9)),s9,'k');xlim([-pi,pi])
xlabel('$t$','interpreter','latex')
ylim([-1 4.7])
set(gca,'FontSize',txtfont,'FontUnits','points')
text(-3,4,'$s_2(t)$','interpreter','latex','FontSize',12)

subplot(223);bar(S4(end/2+2:end/2+5),'k')
xlabel('No. of Harmonic','interpreter','tex')
ylabel('Magnitude','interpreter','tex')
set(gca,'FontSize',txtfont,'FontUnits','points')

subplot(224);bar(S9(end/2+2:end/2+10),'k')
xlabel('No. of Harmonic','interpreter','tex')
ylabel('Magnitude','interpreter','tex')
set(gca,'FontSize',txtfont,'FontUnits','points')

set(gcf,'color','white')

AM = 1 + 0.05*sqrt(t);

px = [0.1307,0.5711];
py = [0.77,0.43,0.1];

MatyLim = [10,10;2,2;2.2,3.3];
Matylim = [-20,-20;-4,-2;0,0];
txtsnr = {'Noiseless','10 dB','0 dB'};
Matpostxt = [7.2,7.2;1.4,1.6;2,2.95];
Xleg = [0.332 0.77];
figure(2)
set(gcf,'Position',[9 197 419 487])
for l=1:length(vD0)
    D0 = vD0(l);
    Nr = 30;
    for m=1:length(vRSR)
        s = 0;
        for i=1:D0
            s = s + amp(i)*cos(2*pi*i*phi);
        end
        s = AM.*s;
        s = s';
        RSR = vRSR(m);
        noi = 10^(-RSR/20)*std(s)*randn(size(s));
        
        s = s + noi;
        
        [F,~] = STFT_Gauss(s,length(s),1e-4,0.5);
        c = ridge_ext(F,0.1,0.1,10,10);
        est_desvGRe= median(abs(real(F(:))))/0.6745;
        est_desvGIm = median(abs(imag(F(:))))/0.6745;
        est_desvG = sqrt(est_desvGRe^2+est_desvGIm^2);
        est_desv = est_desvG/G;
        rmax = floor(0.5*fs/max(c));
        
        vc = [3,5,8,12];
        Wan = zeros(rmax,length(vc));
        
        [F,sF] = STFT_Gauss(s,length(s),1e-4,0.5);
        RSS = zeros(1,rmax);
        GCV = RSS;
        R_l = RSS;
        SCH = RSS;
        Kv = zeros(rmax,H);
        Wn = zeros(rmax,length(vc));
        g = zeros(1,rmax);
        for r=1:rmax
            [s_est,v] = WSF(s,c,r,1,30,F,sF);
            E = (s-s_est);
            RSS(r) = sum(E.^2);
            GCV(r) = N*RSS(r)/(length(s) - 2*r - 1)^2;
            R_l(r) = 1/N*RSS(r) + 2*(est_desv)^2*(2*r + 1)/N;
            for rc=1:length(vc)
                Wn(r,rc) = log10(1/N*RSS(r)) + vc(rc)*r*log10(N)/N;
            end
            for h=1:H
                [~,var] = lpc(E,h);
                Kv(r,h) = log(var) + (5*r + h)*log10(N)/N;
            end
        end
        
        figure(2)
        subplot(3,2,2*(m-1)+l)
        
        [~,bg] = min(GCV);
        plot(GCV,'Linewidth',2)
        hold on
        
        [~,brl] = min(R_l);
        plot(R_l-1,'Linewidth',2)
        
        [~,aux] = min(Wn(:));
        [br,bc] =  ind2sub(size(Wn),aux);
        plot(Wn(:,bc),'Linewidth',2)
        opt_c(l,m) = bc;
        
        [~,aux] = min(Kv(:));
        [br2,bc2] =  ind2sub(size(Kv),aux);
        plot(Kv(:,bc2),'Linewidth',2)
        opt_h(l,m) = bc2;
        plot(bg,GCV(bg),'^r','MarkerFaceColor','r','MarkerSize',5)
        plot(brl,R_l(brl)-1,'^r','MarkerFaceColor','r','MarkerSize',5)
        plot(br,Wn(br,bc),'^r','MarkerFaceColor','r','MarkerSize',5)
        plot(br2,Kv(br2,bc2),'^r','MarkerFaceColor','r','MarkerSize',5) 
        
        
        xlim([0,30])
        set(gca,'XTick',[0,5,10,15,20,25,30])
        ylim([Matylim(m,l),MatyLim(m,l)])
        xlbc = xlabel('Model order','interpreter','tex');
        xlbc.Position(2) = xlbc.Position(2) + 0.1;
        set(gca,'FontSize',txtfont)
        hold off
        xlbc = xlabel('Model order','interpreter','tex');
        xlbc.Position(2) = xlbc.Position(2) + 0.1;
        set(gca,'FontSize',txtfont)
        set(gca,'Position',[px(l),py(m),0.3339,0.22])
        xleg=legend('G','R','W','K','Orientation','vertical','Location','northeast');
        xleg.Position(1) = Xleg(l);
        text(0.5,Matpostxt(m,l),txtsnr{m},'interpreter','tex','FontSize',12,'FontUnits','points')
        
    end
end
[~,b] = min(Wan(:));
[br,~] = ind2sub(size(Wan),b);