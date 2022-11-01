clear all
clc

%% --- Input parameters --------------------------------------------------

t = 0;                      % initial time
T = 1;                      % terminal time
N = 253;                    % daily swap
tj = linspace(t,T,N);       % discretizing time
dj = tj(2)-tj(1);           % stepping time

x0 = 2;                     % initial value x at time t = 0
S0 = -10:1:10;              % initial value δ at time t = 0

p = 2;                      % degree of moment
Ns = 6001;                  % Number of step for MC       
Np = 1000;                  % Number of path for MC

r = 0.05;
kap = 1.187;
sig1 = 0.212;
sig2 = 0.187;
alp = 0.082;
lam = 0.093;
rho = 0.845;

%% --- Our analytical formula --------------------------------------------

tic
for l = 1:length(S0)

    s0 = S0(l);             % vary each initial δ

    Ba = @(tau) (exp(-kap*tau)-1)/kap;
    Ca = @(tau) -(-1+exp(-kap*tau)+kap*tau)*(kap*alp-lam)/kap^2 + tau*(r-1/2*sig1^2);
    Bb = @(tau) exp(-kap*tau);
    Cb = @(tau) (kap*alp-lam)*(1-exp(-kap*tau))/kap;
    Caa = @(tau) tau*sig1^2 - 2*(-1+exp(-kap*tau)+kap*tau)*rho*sig1*sig2/kap^2 + exp(-2*kap*tau)*(-1+4*exp(kap*tau)+exp(2*kap*tau)*(-3+2*kap*tau))*sig2^2/(2*kap^3);
    Cbb = @(tau) 1/2*sig2^2*(1-exp(-2*kap*tau))/kap;
    Qa  = @(tau) x0+Ba(tau)*s0+Ca(tau);
    Qaa = @(tau) Caa;
    Qb  = @(tau) Bb(tau)*s0+Cb(tau);
    Qbb = @(tau) Cbb(tau);

    for k = 1:N-1
        S(k) = 0;
        for j = 0:p
            s1 = 0;
            s2 = 0;
            for i = 0:floor(j/2)
                s1 = s1 + DFactorial(2*i-1)*nchoosek(j,2*i)*Ca(dj)^(j-2*i)*Caa(dj)^i;
            end
            for i = 0:floor((p-j)/2)
                s2 = s2 + DFactorial(2*i-1)*nchoosek(p-j,2*i)*Qb(tj(k)-t)^(p-j-2*i)*Qbb(tj(k)-t)^i;
            end
            S(k) = S(k) + nchoosek(p,j)*Ba(dj)^(p-j)*s1*s2;
        end
    end
    Uex(l) = sum(S)*100^2/T;
end
toc

%% --- Monte Carlo simulation -------------------------------------------

tic
for l = 1:length(S0)
    s0 = S0(l);                 % vary each initial δ
    tt = linspace(t,T,Ns);
    dt = tt(2)-tt(1);
    x(:,1) = x0*ones(Np,1);
    s(:,1) = s0*ones(Np,1);
    for i = 1:Ns-1
        dW1 = sqrt(dt)*randn(Np,1);
        dW2 = sqrt(dt)*randn(Np,1);
        x(:,i+1) = x(:,i) + (r-s(:,i)-1/2*sig1^2)*dt + sig1*dW1;
        s(:,i+1) = s(:,i) + (kap*alp-kap*s(:,i)-lam)*dt + rho*sig2*dW1 + sig2*sqrt(1-rho^2)*dW2;
    end
    M = (Ns-1)/(N-1);
    Umc(l) = mean(sum((x(:,M+1:M:Ns)-x(:,1:M:Ns-M)).^p'))*100^2/T;
end
toc

%% --- Plot graph and display solution -----------------------------------

plot(S0,Uex,'LineWidth',1.5)
hold on
plot(S0,Umc,'o','LineWidth',1.2)

grid on
grid minor

xlabel('$\delta_{t_0}$','Interpreter','latex')
ylabel('$\mathbf{E}_{t_0}^{\mathrm{Q}}\left[\sigma_M^5\right]$','Interpreter','latex')
legend({'$K_M$ with $p=2$ from our formula','$K_M$ with $p=2$ from MC simulations'},'Interpreter','latex','Location', 'northeast')
hold off

[Uex' Umc' abs((Uex'-Umc')./Uex')]          % display solution
mean(abs((Uex'-Umc')./Uex'))                % display error