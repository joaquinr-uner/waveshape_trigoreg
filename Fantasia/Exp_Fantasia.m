%% File that runs the experiment presented in Sec. 6 of of the paper
% Ruiz, J. Colominas, MA. "Waveshape Model Order Estimation by
% Trigonometric Regression". 
% Files from the Fantasia Database are required.

st = 1e5;
ed = 1.5e5-1;
fsi = 250;

N = 1e7;
sigma = 5e-6;
fmax = 0.2;
redun = 1/10;

drt = 'Fantasia_Data';
listing = dir(drt);

for i=3:length(listing)
    fname = listing(i).name;
    data = load([drt '\' fname]);
    s = data.val(2,st:ed)';

    t1 = 0:1/fsi:length(s)/fsi-2/fsi;
    fs = fsi;
    t = 0:1/fs:N/fs-1/fs;
    
    s = s - mean(s);
    N = length(s);
    f = 0:fs/(N*redun):fs*fmax-fs/(N*redun);
    
    fprintf(['Processing Register ' fname '\n'])
    tic;z|
    [F,sF] = STFT_Gauss(s,redun*N,sigma,fmax);
    telp = toc;
    fprintf(['STFT Time Elapsed: ' num2str(telp) '\n'])
    
    tic;
    U = istct_fast(F,f);
    T = toc;
    fprintf(['elapsed time STCT: ' num2str(T) '\n'])
    W = U.*F;
    c = ridge_ext(W,1,1,3,10,redun);
    rmax = floor(fs*0.5/max(f(c)));
    H = floor((log10(N))^2);
    
    RSS = zeros(1,rmax);
    GCV = RSS;
    Rl = RSS;
    Wang = RSS;
    Kav = zeros(rmax,H);
    
    G = sqrt(sqrt(pi)/sqrt(2*sigma));
    est_desvGRe= median(abs(real(F(:))))/0.6745;
    est_desvGIm = median(abs(imag(F(:))))/0.6745;
    est_desvG = sqrt(est_desvGRe^2+est_desvGIm^2);
    est_desv = est_desvG/G;
    b=10;
    for r=1:rmax
        s_est = WSF(s,c,r,1,b,F,sF);
        SE = (s-s_est);
        RSS(r) = sum(SE.^2);
        GCV(r) = N*RSS(r)/(N - 2*r-1)^2;
        Rl(r) = 1/N*RSS(r) + 2*(est_desv)^2*(2*r+1)/N;
        Wang(r) = log10(1/N*RSS(r)) + 2.1*r*log(N)/N;
        for h=1:H
            [~,var] = lpc(SE,h);
            Kav(r,h) = log(var) + (5*r + h)*log(N)/N;
        end
    end
    
    re = findmin(GCV);
    fprintf(['WSF Order ECG: ' num2str(re) '\n'])
    [s_est,v] = WSF(s,c,re,1,b,F,sF);
    
    wsf = cosenos(ones(1,N),linspace(-1/2,1/2,N),re)*v;
    
    resp = data.val(1,st:ed)';
    resp = resp - mean(resp);
    fmaxr = 0.05;
    redunr = 1;
    sigmar = 1e-4;
    [Fresp,sFresp] = STFT_Gauss(resp,redun*N,sigma,fmax);
    
    c_resp = ridge_ext(Fresp,0.1,0.1,10,10);
    c_resp = c_resp+1;
    rmax_resp = floor(fs*0.5/max(f(c_resp)));
    H = floor((log10(N))^2);
    
    RSS = zeros(1,rmax_resp);
    GCV = RSS;
    Rl = RSS;
    Wang = RSS;
    Kav = zeros(rmax_resp,H);
    
    G = sqrt(sqrt(pi)/sqrt(2*sigma));
    est_desvGRe= median(abs(real(Fresp(:))))/0.6745;
    est_desvGIm = median(abs(imag(Fresp(:))))/0.6745;
    est_desvG = sqrt(est_desvGRe^2+est_desvGIm^2);
    est_desv = est_desvG/G;
    
    for r=1:rmax_resp
        resp_est = WSF(resp,c_resp,r,1,5,Fresp,sFresp);
        SE = (resp-resp_est);
        RSS(r) = sum(SE.^2);
        GCV(r) = N*RSS(r)/(N - 2*r-1)^2;
        Rl(r) = 1/N*RSS(r) + 2*(est_desv)^2*(2*r+1)/N;
        Wang(r) = log10(1/N*RSS(r)) + 2.1*r*log(N)/N;
        for h=1:H
            [~,var] = lpc(SE,h);
            Kav(r,h) = log(var) + (5*r + h)*log(N)/N;
        end
    end
    
    re_resp = findmin(GCV);
    
    fprintf(['WSF Order Resp: ' num2str(re_resp) '\n'])
    
    [resp_est,v_resp] = WSF(resp,c_resp,re_resp,1,5,Fresp,sFresp);
    
    wsf_resp = cosenos(ones(1,N),linspace(-1/2,1/2,N),re_resp)*v_resp;
    
    R = struct('ECG_Estimate',s_est,'ECG_Coeff',v,'ECG_WSF',wsf,...
        'Resp_Estimate',resp_est,'Resp_Coef',v_resp,'Resp_WSF',wsf_resp,'rmax'...
        ,rmax,'rmax_resp',rmax_resp);
end
