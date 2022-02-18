load('wave_indexes.mat')

indx = Meta.wave_indexes;
DB ={	
     'AbdAorta';
};

N = 7000;
T = 120;

t = linspace(0,T,N);
SNR = 0:10:30;
NR = length(SNR);
Nr = 70;
SNRout = zeros(Nr,11,NR);
br = zeros(Nr,9,NR);
MMAE = zeros(Nr,9,NR);
Cor = zeros(Nr,9,NR);
CS_est = cell(Nr,length(SNR));
CV = cell(Nr,length(SNR));
fs = N/T;
f = 0:fs/N:(fs/2-fs/N);

sdihr = 0.035;
lfhf = 3;
inf = 0.05;
inf2 = 0.2;

lfhfa = 1/lfhf^2;

phi = lfhf*sdihr/(2*inf*pi)*cos(2*inf*pi*t) + (1-lfhfa)*sdihr/(2*inf2*pi)*cos(2*pi*inf2*t) + 1*t;

for l=1:size(DB,1)
    
    Data = readtable(['PWs_' DB{l} '_PPG.csv'], 'HeaderLines',1);
    
    Data = table2array(Data);
    Data(:,1:2) = [];
    [M,~] = size(Data);
    for j=1:Nr
        
        pw = Data(indx(j),2:end);
        pw = pw';
        pw(isnan(pw)) = [];
        pw = pw - mean(pw);
        pw = pw/max(pw);
        sh = fft(pw);
        sh = sh(2:floor(end/2));
        Nf = length(sh);
        
        v = [2*real(sh); -2*imag(sh)];
        sigma = 3.33e-4;
        sr = cosenos(1+0.02*sqrt(t),phi,Nf)*v;
        sr = sr - mean(sr);
        src = sr;
        redun = 1;
        for k=1:length(SNR)
            sr = src + 10^(-SNR(k)/20)*std(src)*randn(size(src));
            
            H = 1;
            [F,sF] = STFT_Gauss(sr,length(sr)*redun,sigma,0.5);
            c = ridge_ext(F,0.1,0.1,50,10,redun);
            rmax = floor(length(sr)*0.5/max(c));
            RSS = zeros(1,rmax);
            GCV = RSS;
            Rl = RSS;
            SCH = RSS;
            Wang = RSS;
            Kav = zeros(rmax,H);
            G = sqrt(sqrt(pi)/sqrt(2*sigma));
            est_desvGRe= median(abs(real(F(:))))/0.6745;
            est_desvGIm = median(abs(imag(F(:))))/0.6745;
            est_desvG = sqrt(est_desvGRe^2+est_desvGIm^2);
            est_desv = est_desvG/G;
            for r=1:rmax
                s_est = WSF(sr,c,r,1,50,F,sF);
                SE = (sr-s_est);
                RSS(r) = sum(SE.^2);
                GCV(r) = N*RSS(r)/(N - 2*r)^2;
                Rl(r) = 1/N*RSS(r) + 2*(est_desv)^2*(2*r)/N-est_desv;
                Wang(r) = log10(1/N*RSS(r)) + 2.1*r*log(N)/N;
                for h=1:H
                    [~,var] = lpc(SE,h);
                    Kav(r,h) = log(var) + (5*r + h)*log(N)/N;
                end
            end
                      
            [~,br(j,1,k)] = min(GCV);
            [~,br(j,2,k)] = min(Rl);
            [~,br(j,4,k)] = min(Wang);
            [~,b] = min(Kav(:));
            [br(j,5,k),~] = ind2sub(size(Kav),b);
            br(j,6:9,k) = [3 6 9 18];
            S_est = zeros(length(sr),size(br,2)+2);
            V_est = cell(size(br,2));
            for i=1:size(br,2)
                [S_est(:,i),V_est{i}] = WSF(sr,c,br(j,i,k),1,50,F,sF);
                SNRout(j,i,k) = 20*log10(std(src)/std(src-S_est(:,i)));
                pw_est = cosenos(ones(size(pw)),linspace(-1/2,1/2,length(pw))',br(j,i,k))*V_est{i};
                pw_est = pw_est/max(pw_est);
                
                pwcorr = xcorr(pw_est,pw,length(pw));
                autcorrpw = xcorr(pw,pw,0);
                [mpwcorr,tau] = max(pwcorr);
                tau = length(pw) - tau;
                pw_est2 = circshift(pw_est',[tau,0]);
                
                cor = round(dot(pw_est2,pw)/dot(pw,pw)*100);
                Cor(j,i,k) = cor;
                
                MAE = sum(abs(pw_est2'-pw))/sum(abs(pw));
                MMAE(j,i,k) = MAE;
            end
            
            F_hthr = threshold(F,1);
            F_sthr = threshold(F,2);
            S_est(:,10)= 2/max(sF)*real(sum(F_hthr,1));
            SNRout(j,10,k) = 20*log10(std(src)/std(src-S_est(:,10)));
            
            S_est(:,11) = 2/max(sF)*real(sum(F_sthr,1));
            SNRout(j,11,k) = 20*log10(std(src)/std(src-S_est(:,11)));
            
            CS_est{j,k} = S_est;
            CV{j,k} = V_est;
        end
    end
    Meta = struct('wave_indexes',indx,'SNR',SNR,'fixed_r',br,'Correlation',Cor,'MAE',MMAE,'SNRout',SNRout,...
        'Estimations',CS_est,'Coeffs',CV,'waveindx',indx);
    save(['Results_' DB{l} 'WithThr-17-07-20.mat'],'Meta')
end
