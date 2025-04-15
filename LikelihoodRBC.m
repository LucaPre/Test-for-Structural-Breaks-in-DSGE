% Computes the (negative) log Likelihood of the full sample with the Kalman filter
% assuming a switch in parameters at t=breakdate+1 (Technical Note: constant term of the log Likelihood due to
% integrating constant is left out)
function [crit, states] = LikelihoodRBC(theta,data,breakdate)
% parameters that are used for optimization
omega1=theta(1);
rho1=theta(2);

omega2=theta(3);
rho2=theta(4);

% remaining parameters that are correctly calibrated
alpha=0.3;
beta=0.99;
delta=0.02;

kstar=((1/beta-1+delta)/alpha)^(1/(alpha-1)); % steady state k
cstar=kstar^alpha+(1-delta)*kstar-kstar; % steady state c

deltaderiv=0.0001; % small dx for numerical derivative
[Hy1, Hx1] = SolveFOP(alpha,beta,delta,rho1,omega1,deltaderiv); % solve the model before break
[Hy2, Hx2] = SolveFOP(alpha,beta,delta,rho2,omega2,deltaderiv); % solve the model after break

% State space and Kalman filter notation from 2023 Continuing Education Webcasts: Techniques of Empirical Macroeconomics—Òscar Jordà and Karel Mertens, slides part 2 
% st=[1 kt zt]'
F=[1 0 0; kstar-Hx1(1)*kstar Hx1(1) Hx1(2);0 0 rho1];
Qroot=[0 0 0; 0 0 0; 0 0 omega1];
Q=Qroot*Qroot';
H=[cstar-Hy1(1)*kstar Hy1(1) Hy1(2)];
R=0;

state_filt=[1;kstar;0]; % Initial state
P_stat_vec=(eye(4)-kron(F(2:end,2:end),F(2:end,2:end)))^-1*[zeros(3,1);Q(3,3)]; % Initial variance from stationary distribution
P_filt(2:3,2:3)=reshape(P_stat_vec,2,2);
states=zeros(length(data),3);
logL=0;

% Recursion until breakdate
for t=1:breakdate
state_fore=F*state_filt; % State Forecast mean
P_fore=F*P_filt*F'+Q; % State Forecast Variance
y_fore=H*state_fore; % Observation Forecast mean
G_fore=H*P_fore*H'+R; % Observation Forecast Variance
state_filt=state_fore+P_fore*H'*G_fore^-1*(data(t)-y_fore); % filtered State estimate
states(t,:)=state_filt';
P_filt=P_fore-P_fore*H'*G_fore^-1*H*P_fore; % Variance of filtered State
logL=logL-0.5*log(det(G_fore))-0.5*(data(t)-y_fore)'*G_fore^-1*(data(t)-y_fore); % log likelihood summation (full sample likelihood equivalent to product of likelihoods of yt conditional on past observations)
end

% Update state space matrices after break
F=[1 0 0; kstar-Hx2(1)*kstar Hx2(1) Hx2(2);0 0 rho2];
Qroot=[0 0 0; 0 0 0; 0 0 omega2];
Q=Qroot*Qroot';
H=[cstar-Hy2(1)*kstar Hy2(1) Hy2(2)];

% Recursion after breakdate
for t=breakdate+1:length(data)
state_fore=F*state_filt; % State Forecast mean
P_fore=F*P_filt*F'+Q; % State Forecast Variance
y_fore=H*state_fore; % Observation Forecast mean
G_fore=H*P_fore*H'+R; % Observation Forecast Variance
state_filt=state_fore+P_fore*H'*G_fore^-1*(data(t)-y_fore); % filtered State estimate
states(t,:)=state_filt';
P_filt=P_fore-P_fore*H'*G_fore^-1*H*P_fore; % Variance of filtered State
logL=logL-0.5*log(det(G_fore))-0.5*(data(t)-y_fore)'*G_fore^-1*(data(t)-y_fore); % log likelihood summation (full sample likelihood equivalent to product of likelihoods of yt conditional on past observations)
end

crit=-logL; % Return criterion (negative of likelihood because we're using minimization algorithms)

% Introduce bounds by penalty
if omega1<0
crit=1e+50;
elseif rho1<0
crit=1e+50;
elseif rho1>1
crit=1e+50;
elseif omega2<0
crit=1e+50;
elseif rho2<0
crit=1e+50;
elseif rho2>1
crit=1e+50;
end
