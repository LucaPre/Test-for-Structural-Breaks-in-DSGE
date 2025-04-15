function f = scores(H,deltaderiv,alpha,beta,delta,rho,omega)
kstar=((1/beta-1+delta)/alpha)^(1/(alpha-1)); % steady state k
Hy=H(1,:);
Hx=H(2,:);

[f11, f21]=egamma(kstar+deltaderiv,0,0,alpha,beta,delta,rho,omega,Hy,Hx);
[f12, f22]=egamma(kstar-deltaderiv,0,0,alpha,beta,delta,rho,omega,Hy,Hx);
f(1:2)=([f11 f21]-[f12 f22])/(2*deltaderiv); % Derivative of gamma vector wrt to kt around steady state

[f11, f21]=egamma(kstar,0+deltaderiv,0,alpha,beta,delta,rho,omega,Hy,Hx);
[f12, f22]=egamma(kstar,0-deltaderiv,0,alpha,beta,delta,rho,omega,Hy,Hx);
f(3:4)=([f11 f21]-[f12 f22])/(2*deltaderiv); % Derivative of gamma vector wrt to zt around steady state

[f11, f21]=egamma(kstar,0,0+deltaderiv,alpha,beta,delta,rho,omega,Hy,Hx);
[f12, f22]=egamma(kstar,0,0-deltaderiv,alpha,beta,delta,rho,omega,Hy,Hx);
f(5:6)=([f11 f21]-[f12 f22])/(2*deltaderiv); % Derivative of gamma vector wrt to sigma around steady state




