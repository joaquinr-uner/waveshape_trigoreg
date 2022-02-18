addpath(genpath('/home/jruiz'))

res = 600;
D0 = [4,3];
vRSR = [0,10,20];
N = 3000;

t = linspace(0,1,N);
cw = 2.1;

fs = 1/(t(2)-t(1));

G = sqrt(sqrt(pi)/sqrt(2*1e-4));

a1 = 120;
b1 = 15;
phi1 = a1*t  + b1*t.^2; Trig = 'Pure';
f1 = a1 + 2*b1*t;

f02 = 200;
fa = 25;
%phi2 = a2*t + b2*t.^2;
phi2 = f02*t + fa/(2*pi)*cos(2*pi*t);
f2 = f02 - fa*sin(2*pi*t);

%A1 = 1 + 0.05*sqrt(t);
A1 = ones(1,N);
%A2 = 1 + 0.1*sin(2*pi*t);
A2 = ones(1,N);

Nr = 1000;

for i=1:length(vRSR)
    DG = zeros(2,Nr);
    DR = DG;
    DW = DG;
    DK = DG;

    Mamp = {zeros(D0(1),Nr),zeros(D0(2),Nr)};
    for j=1:Nr
        amp1 = [1; 0.8*rand(D0(1)-1,1)+0.1];
        Mamp{1}(:,j,i) = amp1;
        amp2 = [0.95; 0.8*rand(D0(2)-1,1)+0.1];
        Mamp{2}(:,j,i) = amp2;
        
        s1 = Fdict(A1,phi1,length(amp1))*[amp1; zeros(size(amp1))];
        s2 = Fdict(A2,phi2,length(amp2))*[amp2; zeros(size(amp2))];
        sc = s1 + s2;
        n = 10^(-vRSR(i)/20)*std(sc)*randn(length(sc),1);
        
        s = sc + n;
        
        [F, sF] = STFT_Gauss(s, length(s),1e-4,0.5);
        
        est_desvGRe= median(abs(real(F(:))))/0.6745;
        est_desvGIm = median(abs(imag(F(:))))/0.6745;
        est_desvG = sqrt(est_desvGRe^2+est_desvGIm^2);
        est_desv = est_desvG/G;
        
        c1 = ridge_ext(F,0.1,0.1,10,10);
        
        y1 = zeros(N,1);
        b = 20;
        for k=1:length(y1)
            binf = max([1,c1(k)-b]);
            bup = min([size(F,1),c1(k)+b]);
            y1(k) = 1/max(sF)*sum(F(binf:bup,k));
        end
        
        s2_est = s - 2*real(y1);
        
        
        [F2_est, sF2] = STFT_Gauss(s2_est,length(s2_est),1e-4,0.5);
        
        c2 = ridge_ext(F2_est,0.1,0.1,10,10);
        
        if mean(c1>c2)
            caux = c1;
            c1 = c2;
            c2 = caux;
        end
        
        rmax1 = floor(N/2/max(c1));
        rmax2 = floor(N/2/max(c2));

        RSS = zeros(rmax1,rmax2);
        GCV = RSS;
        Wn = RSS;
        R_l = RSS;
        Kv = RSS;
        h = floor(log(N)^2);
        for r1=1:rmax1
            for r2=1:rmax2
                s_est = WSF2(s,c1,c2,r1,r2,1,b,F,sF);
                E = (s-s_est);
                RSS(r1,r2) = sum(E.^2);
                GCV(r1,r2) = N*RSS(r1,r2)/(length(s) - 2*(r1+r2) - 1)^2;
                Wn(r1,r2) = log(1/N*RSS(r1,r2)) + cw*(r1+r2)*log(N)/N;
                R_l(r1,r2) = 1/N*RSS(r1,r2) + 2*(est_desv)^2*(2*(r1+r2) + 1)/N;
                [~,varK] = lpc(E,h);
                Kv(r1,r2) = log(varK) + (5*(r1+r2) + h)*log(N)/N;
            end
        end
        
        [~,b1] = min(GCV(:));
        [rg1, rg2] = ind2sub(size(GCV),b1);
        DG(:,j) = [rg1,rg2];
        
        [~,b2] = min(R_l(:));
        [rr1,rr2] = ind2sub(size(R_l),b2);
        DR(:,j) = [rr1,rr2];

        [~,b3] = min(Wn(:));
        [rw1,rw2] = ind2sub(size(Wn),b3);
        DW(:,j) = [rw1,rw2];
        
        [~,b4] = min(Kv(:));
        [rk1,rk2] = ind2sub(size(Kv),b4);
        DK(:,j) = [rk1,rk2];
        fprintf(['Iteracion: ' num2str(j) ' de ' num2str(Nr) ' RSR: ' num2str(vRSR(i)) '\n'])
    end
    
    
    SA = struct('D0', D0, 'DGCV', DG, 'DRl', DR,'DWang', DW, 'DKavalieris', DK,...
        'matrix_phi', [phi1;phi2], 'matrix_A',[A1;A2],...
        'ampl_coefs1',Mamp{1},'ampl_coefs2',Mamp{2},'n_realiz',Nr);
    
    save(['Results_ExpTrigo_2Chirps_' date '_RSR_' num2str(vRSR(i)) 'dB_LinSin.mat'],'SA')
    
end
