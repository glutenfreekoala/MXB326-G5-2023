% initialise workspace:
clear; close all;

fprintf("Calculating semi-analytical model.\n")

% set constant model parameters
Swr = 0.0375;
Sor = 0.15;
F = 2;
beta_list = [1, 2, 4, 10];

X = linspace(0,1,10);
T = [0, 0.01, 0.05, 0.1, 0.3, 0.8125, 1];

N = length(X);
M = length(T);

S_SemiAnalytical = zeros(N,M,length(beta_list));

% semi analytical solution:
tiledlayout(2,2);
for beta_choice = 1:length(beta_list)
    % semi-analytical solution
    % set domain intervals and discretisation points:
    
    S = zeros(N,M);


    % set model parameters
    beta = beta_list(beta_choice);
    fprintf("Beta = %d\n", beta)
    [gamma, alpha, omega] = paramsFunc(F,beta,Sor,Swr);
     
    % determine x tilde partial derivative of phi
    syms xt t
    phi_dot = diff(phiFunc(xt, t, beta, omega), xt);
    phi_xt = @(xt, t) eval(phi_dot);
    
    % initial condition
    S(:,1) = 1-Swr;
    
    
    % compute S at each x, t=t2:tM
    for n = 2:length(T)
        for i = 1:length(X)
            xtilde = (beta*S(i,n-1) + gamma) * X(i) + omega*T(n-1); % guess
    
            xtilde = fzero(@(x) phiFunc(x, T(n), beta, omega) - exp(alpha*X(i)), xtilde); % find xtilde
            S(i, n) = (1/beta)*((alpha * exp(alpha * X(i)))/phi_xt(xtilde, T(n))-gamma); % calc S at xi tn
        end
        fprintf("Completed calculation of timestep %d \n", n)
    end

    S_SemiAnalytical(:,:,beta_choice) = S;

    nexttile; hold on; % plot S(x,t) for selected beta at each time
    for n=2:7
        plot(X,S(:,n))
    end
    xlim([0 1]); ylim([0 1]);
    title("Beta = "+num2str(beta_list(beta_choice))); ylabel("Oil Saturation"); xlabel("x");
end

% save values for semi analytical solution and clean workplace.
X_SemiAnalytical = X;
T_SemiAnalytical = T;
clear xt t X S N M xtilde beta beta_choice i n phi_dot phi_xt ;

% Numerical solution:



% Analysis:


% Average Saturation


% Buckley-Leverett
fprintf("Plot semi-analytical solution with buckley leverett solution")
figure; hold on;
X = linspace(0,1,1000);
S_BuckleyLeverett = BuckleyLeverett(X, T_SemiAnalytical, Swr, Sor);
plot(X, S_BuckleyLeverett, "--"); plot(X_SemiAnalytical, S_SemiAnalytical(:,:,4));
xlim([0 1]); ylim([0 1]); xlabel("x"); ylabel("Saturation");
title("Semi Analytical Solution with Buckley-Leverett Shock")


