function [Hy, Hx] = SolveFOP(alpha,beta,delta,rho,omega,deltaderiv)

% Starting values close to actual values with rho=0.95,0.85 for faster computation with 
% benefit of hindsight, stability conditions in this coding exercise are 
% always fulfilled but should generally be checked additionally as there exist 2 solutions

Hystart=[0.05 0.5 0]; % Start value for Hy 
Hxstart=[0.95 2.15 0]; % Start value for Hx
eigval=2;
while eigval>1 % Check stability condition
H=[Hystart;Hxstart];
fixedFunction=@(x) scores(x,deltaderiv,alpha,beta,delta,rho,omega);
Hsolved=fsolve(fixedFunction,H); % Solve for Hy and Hx so that derivatives wrt to kt, zt and sigma are 0
Hy=Hsolved(1,:);
Hx=Hsolved(2,:);
eigval=abs(Hx(1)); % Largest Eigenvalue (rho is already restricted between 0 and 1 by likelihood function)

% Update starting values randomly if stability condition not satisfied
Hystart=Hystart+[mvnrnd(zeros(1,2),eye(2)) 0]; 
Hxstart=Hxstart+[mvnrnd(zeros(1,2),eye(2)) 0];
if eigval>1
    disp("Warning: Solution unstable")
end
end

