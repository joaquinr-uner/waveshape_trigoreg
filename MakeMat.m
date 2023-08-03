function [M,L,marg1,marg2] = MakeMat(D)
    [K,N] = size(D);

    L = zeros(K,2);
    for k=1:K
        L(k,1) = min(D(k,:));
        L(k,2) = max(D(k,:));
    end
    
    M = zeros(L(1,2),L(2,2));
    for n=1:N
       n1 = D(1,n);
       n2 = D(2,n);
       M(n1,n2) = M(n1,n2) + 1;
    end
    marg1 = hist(D(1,:),1:L(1,2));
    marg2 = hist(D(2,:),1:L(2,2));
end


