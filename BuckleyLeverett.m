function S = BuckleyLeverett(x,t,Swr,Sor)

v = 1/(1-Swr-Sor);
S=zeros(length(x),length(t));

for n = 1:length(t)
    for i = 1:length(x)
        if x(i) < v*t(n)
            S(i,n) = Sor;
        else
            S(i,n) = 1-Swr;
        end
    end
end
end