% Buckley-Leverett Shock Solution

X = linspace(0,1,1000);
T = [0, 0.01, 0.05, 0.1, 0.3, 0.8125, 1];

Swr = 0.0375;
Sor = 0.15;

S_Shock = BuckleyLeverett(X,T,Swr,Sor);

figure; hold on;
ylabel("Oil Saturation"); xlabel("x"); title("Buckley-Leverett Approximate solution for large beta");
for n = 1:length(T)
    plot(X, S_Shock(:,n))
end
    