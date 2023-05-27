% initialise workspace
clear; close all;
tiledlayout(2,2);

beta_list = [1, 2, 4, 10];
for beta_choice = 1:length(beta_list)

    % semi-analytical solution
    % set domain intervals and discretisation points:
    X = linspace(0,1,1000);
    T = [0, 0.01, 0.05, 0.1, 0.3, 0.8125, 1];
    
    N = length(X);
    M = length(T);
    
    S = zeros(N,M+1);


    % set model parameters
    beta = beta_list(beta_choice);
    Swr = 0.0375;
    Sor = 0.15;
    F = 2;
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
    end

    nexttile; hold on; % plot S(x,t) for selected beta at each time
    plot(X,S(:,2))
    plot(X,S(:,3))
    plot(X,S(:,4))
    plot(X,S(:,5))
    plot(X,S(:,6))
    plot(X,S(:,7))
    xlim([0 1]); ylim([0 1]);
    title("Beta = "+num2str(beta_list(beta_choice))); ylabel("Oil Saturation"); xlabel("x");
end