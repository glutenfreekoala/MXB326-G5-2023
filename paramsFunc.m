function [gamma, alpha, omega] = paramsFunc(F, beta, Sor, Swr)

% Solve eqn (6) for gamma
gamma = beta/(F-1) * (1-Swr-F*Sor);

% Use gamma to solve eqn (11) for alpha
alpha = -beta^2 * ((Sor + gamma/beta)*(1-Swr+ gamma/beta))/(1-Swr-Sor);

% Use gamma, alpha to solve eqn (12) for omega
omega = beta*alpha / (beta-beta*Swr + gamma);
