% Computes canonical form with plugged in policy functions
function [f1, f2] = egamma(K,Z,Sigma,alpha,beta,delta,rho,omega,Hy,Hx)
kstar=((1/beta-1+delta)/alpha)^(1/(alpha-1)); % steady state k
cstar=kstar^alpha+(1-delta)*kstar-kstar; % steady state c
zstar=0; % steady state z
nut1=0; % steady state shock 
C=cstar+Hy*[K-kstar;Z-zstar;Sigma]; % linear policy function of ct
Kt1=kstar+Hx*[K-kstar;Z-zstar;Sigma]; % linear policy function of kt+1
Zt1=rho*Z+Sigma*omega*nut1; % evolution of state vector 
Ct1=cstar+Hy*[Kt1-kstar;Zt1-zstar;Sigma]; % linear policy function of ct+1
f1=1/C-beta*1/Ct1*((1-delta)+exp(Zt1)*alpha*Kt1^(alpha-1)); % first row of gamma vector
f2=C+Kt1-(1-delta)*K-exp(Z)*K^alpha; % second row of gamma vector

