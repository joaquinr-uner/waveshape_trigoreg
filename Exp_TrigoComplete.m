%% This file runs the experiment described in Sec. 2.2 of the paper
%% Ruiz, J. Colominas, MA. "Waveshape Model Order Estimation by
%% Trigonometric Regression". 
SNR = [0;5;10;15];
M = size(SNR,1);

VD0 = [1,3,4,6,9];
nD0 = length(VD0);

Nr = 1000;

DoutRl = zeros(M,Nr);
DoutGCV = DoutRl;
DoutWang = DoutRl;
DoutKavalieris = DoutRl;
DoutSCH = DoutRl;

t = linspace(0.001,1,2000);

N = length(t); 

fs = 1/(t(2)-t(1));

phi = 100*t;

phi_2 = 50*t + 15*t.^2;

phi_3 = 60*t + 15/(2*pi)*cos(2*pi*t);

PHI = [phi;phi_2;phi_3];

Type = ['Pure';'CLin';'CSin'];
nT = size(Type,1);

G = sqrt(sqrt(pi)/sqrt(2*1e-4));

vc = [1,2,5,8,12];
H = floor(log10(N)^2);

for i1=1:nD0
    D0 = VD0(i1);
    for i2=1:nT
        phi = PHI(i2,:);
        T = Type(i2,:);
        Mamp = zeros(M,D0,Nr);
        for i3=1:M
            RSR = SNR(i3);
            for j=1:Nr
                tic
                amp = zeros(1,D0);
                amp(1) = 1;
                for i=2:D0
                    amp(i) = 0.8*rand(1,1)+0.1;
                end
                Mamp(i3,:,j) = amp;
                
                x_c = 0;
                for i=1:D0
                    x_c = x_c + amp(i)*cos(2*pi*i*phi);
                end
                x_c = x_c';
                desv_n = 10^(-RSR/20)*std(x_c);
                w_opt = 1e-4;
                
                n = desv_n*randn(size(x_c));
                x = x_c + n;
                redun = 1;
                [F,sF] = STFT_Gauss(x,length(x)*redun,w_opt,0.5);
                c = ridge_ext(F,0.1,0.1,10,10,redun);
                l_max = floor(size(F,1)/(max(c)-1));
                l_max = max([l_max 1]);
                
                est_desvGRe= median(abs(real(F(:))))/0.6745;
                est_desvGIm = median(abs(imag(F(:))))/0.6745;
                est_desvG = sqrt(est_desvGRe^2+est_desvGIm^2);
                est_desv = est_desvG/G;
                
                RSS = zeros(1,l_max);
                GCV = zeros(1,l_max);
                R_l = zeros(1,l_max);
                Wn = zeros(l_max,length(vc));
                Kv = zeros(l_max,H);
                for r=1:l_max
                    [x_recon, v] = WSF(x,c,r,1,50,F,sF);
                    E = (x-x_recon);
                    RSS(r) = sum(E.^2);
                    GCV(r) = N*RSS(r)/(length(x) - 2*r - 1)^2;
                    R_l(r) = 1/N*RSS(r) + 2*(est_desv)^2*(2*r + 1)/N;
                    for rc=1:length(vc)
                            Wn(r,rc) = log10(1/N*RSS(r)) + vc(rc)*r*log(N)/N;
                    end
                    for h=1:H
                        [~,var] = lpc(E,h);
                        Kv(r,h) = log(var) + (5*r + h)*log(N)/N;
                    end
                end
                
%                fprintf(['D: ' num2str(D0) ' Case: ' T ' RSR: ' num2str(RSR) 'dB Realization ' num2str(j) ' out of ' num2str(Nr) '\n'])
                [~,b] = min(R_l);
                DoutRl(i3,j) = b;
%                fprintf(['Rl: ' num2str(b)])
                [~,b] = min(GCV);
                DoutGCV(i3,j) = b;
%                fprintf([' GCV: ' num2str(b)])
                [~,b] = min(Wn(:));
                [br,~] = ind2sub(size(Wn),b);
                DoutWang(i3,j) = br;
%                fprintf([' Wang: ' num2str(br)])
                [~,b] = min(Kv(:));
                [br,~] = ind2sub(size(Kv),b);
                DoutKavalieris(i3,j) = br;
                telp = toc;
%                fprintf([' Kavalieris: ' num2str(br) ' out of ' num2str(l_max)...
%                ' possible candidates (' num2str(SNR(i3)) ' dB). Elapsed Time: ' num2str(telp) '\n'])
            end
        end
        aux = struct('SNR', SNR, 'Type', T, 'n_realiz', Nr,...
                        'D0', D0, 'DGCV', DoutGCV, 'DRl', DoutRl,...
                        'DWang', DoutWang, 'DKavalieris',...
                        DoutKavalieris,'inst_phase',...
                        phi, 'ampl_coefs', Mamp);
            
        save(['Exp_' T '_' num2str(D0) '_Harm.mat'],'aux')

    end
end