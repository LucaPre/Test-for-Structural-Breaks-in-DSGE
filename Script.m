clear
clc
%% Recursive form with first order perturbation
alpha=0.3;
beta=0.99;
delta=0.02;
omega=0.01;
deltaderiv=0.0001; % small dx for numerical derivative

% Solve for policy functions 
[Hy1, Hx1] = SolveFOP(alpha,beta,delta,0.95,omega,deltaderiv); % with rho=0.95
[Hy2, Hx2] = SolveFOP(alpha,beta,delta,0.85,omega,deltaderiv); % with rho=0.85

% Note: You may test the calibration from Fernandez-Villaverde et al.
% (2016), page 42 to confirm that this code gives the same results 

%% Monte Carlo
% running this section takes around 8 minutes at my laptop with 4 cores
MC=500; % Number of Monte Carlo simulations
rng(1) 
seeds=randi([1 2^32-1],MC,1); 
reject=zeros(MC,2);
level=0.001; % Significance level
kstar=((1/beta-1+delta)/alpha)^(1/(alpha-1)); % Steady state k
cstar=kstar^alpha+(1-delta)*kstar-kstar; % Steady state c
Tburn=500; % Burn period
rhos=[0.95 0.85]; % Change second value to 0.95 to simulate under the null of no break
Tpreabreaks=[50 100]; % Sample size after burn period before break
Tpostbreaks=[50 100]; % Sample size after break
for j=1:2
Tprebreak=Tpreabreaks(j);
Tpostbreak=Tpostbreaks(j);
parfor m=1:MC
rng(seeds(m))
rho1=rhos(1);
rho2=rhos(2);
T=Tburn+Tprebreak+Tpostbreak;
Z=zeros(T+1,1);K=zeros(T+1,1);C=zeros(T,1);I=zeros(T,1);Y=zeros(T,1);K(1)=kstar;
% Simulation for burn period and before the break
for t=1:Tburn+Tprebreak
Z(t+1)=rho1*Z(t)+omega*randn(1); % Generate Zt+1
K(t+1)=kstar+Hx1(1:2)*[K(t)-kstar; Z(t)]; % Generate Kt+1 from linear solution
C(t)=cstar+Hy1(1:2)*[K(t)-kstar; Z(t)]; % Generate Ct from linear solution
I(t)=K(t+1)-(1-delta)*K(t); % Generate It
Y(t)=C(t)+I(t); % Generate Yt
end

% Simulation after break
for t=Tburn+Tprebreak+1:T
Z(t+1)=rho2*Z(t)+omega*randn(1);
K(t+1)=kstar+Hx2(1:2)*[K(t)-kstar; Z(t)];
C(t)=cstar+Hy2(1:2)*[K(t)-kstar; Z(t)];
I(t)=K(t+1)-(1-delta)*K(t);
Y(t)=C(t)+I(t);

% For the case rho2=rho1 (if one wants to test the size of the test)
if rho2==rho1
K(t+1)=kstar+Hx1(1:2)*[K(t)-kstar; Z(t)];
C(t)=cstar+Hy1(1:2)*[K(t)-kstar; Z(t)];
I(t)=K(t+1)-(1-delta)*K(t);
Y(t)=C(t)+I(t);
end

end
data=C(Tburn+1:T); % Throw away burn period and use Ct as observation variable

% MLE Without break
theta0=[omega 0.9]; % Starting values for optimization (average of both rhos)
fixedFunction=@(x) LikelihoodRBCNoBreak(x,data);
[thetahat, fval, exitflag]=fminsearch(fixedFunction,theta0);
LogLikFull=-fval;

% MLE With break
theta0=[thetahat thetahat]; % Starting values for optimization centered around estimate without break
fixedFunction=@(x) LikelihoodRBC(x,data,Tprebreak);
[thetahat, fval, exitflag]=fminsearch(fixedFunction,theta0); 
LogLikBreak=-fval;

LR=2*(LogLikBreak-LogLikFull); % Likelihood ratio test statistics
if LR>chi2inv(1-level,2)
    reject(m,j)=1; % Reject if statistic larger than critical value 
end
end
end

disp(['Rejection rate with sample size 50 before and after break is ' num2str(mean(reject(:,1))) ' at the ' num2str(level) ' level']);
disp(['Rejection rate with sample size 100 before and after break is ' num2str(mean(reject(:,2))) ' at the ' num2str(level) ' level']);



