% testing out whether the IVS by Dietz can be implemented as a first-order
% filter

% t = 0.001:0.001:1;
% timeind = 1:1000;
clear

dtau = 0.0001;
tau = dtau:dtau:1;

IPD = (pi*ones(1,10000));
IPD(5000:5999) = linspace(pi,0,1000);
IPD(6000:end) = 0;
IPD = IPD + 0.5*(rand(1,10000)-0.5);
% IPD = pi*rand(1,10000);
IPD = linspace(0,10,length(tau)).*sin(20*tau);

tau_s = 0.04;
alpha = exp(-dtau/tau_s);

IVSG_fof = exp(1i*IPD(1));
IVSG_int = IVSG_fof;
for timInd=2:10000
    tauInd = 1:timInd-1;
    ivsgfun = @(tau)exp(1i*IPD(timInd-tauInd)).*exp(-tau(tauInd)/tau_s);
    IVSG_int(timInd) = (1/tau_s) * (simpleintegral(ivsgfun,tau(tauInd),dtau));
    IVSG_fof(timInd) = (alpha*IVSG_fof(timInd-1)+(1-alpha)*exp(1i*IPD(timInd)));
end

figure
hold on
% plot(exp(-tau/tau_s))
% plot(IVSG)
plot(IPD)
plot(abs(IVSG_int))
plot(abs(IVSG_fof))
% plot(angle(IVSG_int)-angle(IVSG_fof))

%% to show how integration works
xmin = 0;
xmax = 5;
dx = 0.001;
fun = @(x)x;

int_riemann = simpleintegral(fun,xmin:dx:xmax,dx);
int_numeric = integral(fun,xmin,xmax);

%% function definitions
function out = simpleintegral(fun,x,dx)
    out = sum(fun(x))*dx;
end

