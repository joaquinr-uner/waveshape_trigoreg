function [c] = ridge_ext (F,a,b,e,P,redun)
% Function that implements the T-F representation ridge extraction
% procedure discussed in Meignen, T. Oberlin, S. McLaughlin, 
% "A new algorithm for multicomponent signals analysis based on 
% synchrosqueezing: With an applicationto signal sampling and denoising"


if nargin<6
   redun = 1; 
end
[K, N] = size(F);

ti = randi([P,round(N/P)]); %indice de inicio

v_t = round(linspace(ti,N-ti,P));

C = zeros(P,N);
Fun = zeros(P,1);

for p=1:P
   [~,C(p,v_t(p))] = max(abs(F(:,v_t(p))).^2);
   [~,C(p,v_t(p)-1)] = max(abs(F(:,v_t(p)-1)).^2);
   
   for i=v_t(p)+1:N
   I = C(p,i-1)-e:C(p,i-1)+e;
   I(I<1) = 1;
   I(I>K) = K;
   [Fun_aux,aux] = max(abs(F(I,i)).^2 - a*(I'-C(p,i-1)).^2 - b*(I' - 2*C(p,i-1) + C(p,i-2)).^2);
   C(p,i) = aux + C(p,i-1) - e;
   if (C(p,i)<1)
       C(p,i)=1;
   end
   if (C(p,i)>K)
               C(p,i) = K;
   end
   Fun(p) = Fun(p) + Fun_aux;
   end
   
   for i=v_t(p)-1:-1:1
       I = C(p,i+1)-e:C(p,i+1)+e;
       I(I<1) = 1;
       I(I>K) = K;
       [~,aux] = max(abs(F(I,i)).^2 - a*(I'-C(p,i+1)).^2 - b*(I' - 2*C(p,i+1) + C(p,i+2)).^2);
       C(p,i) = aux + C(p,i+1) - e;
       if (C(p,i)<1)
           C(p,i)=1;
       end
       if (C(p,i)>K)
               C(p,i) = K;
       end
       Fun(p) = Fun(p) + Fun_aux;
   end
end

[~,indx] = max(Fun);
c = C(indx,:);
c = c - 2;
end