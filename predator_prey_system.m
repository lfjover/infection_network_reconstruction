
%% equation parameters
clear all close all
nH = 20; %number of Plants
nV = 20; %number of Animals
tfinal = 200;
dt = 0.01;
t = dt:dt:tfinal;

% load('matrices/matrices_10by10_100_invertible.mat')
% M = matrices(:,:,100);
 M = double(rand(nH,nV)>0.98);
 M = M - diag(diag(M)) + eye(size(M));


rMin = 0.15; rMax=1.8;
phiMin = 10^-8; phiMax = 10^-7;
betaMin = 10; betaMax = 20;
mMin = 0.0015; mMax = 0.02;
K = mMax/betaMin/phiMin*10;
Hmin = 10^4; Hmax = 10^5;
Vmin = 10^5; Vmax = 10^6;
a = ones(nH);

phiStar = phiMin + (phiMax - phiMin)*rand(1,nV);
betaStar = betaMin + (betaMax - betaMin)*rand(1,nV);
betaStarM = repmat(betaStar,nH,1);
phiStarM = repmat(phiStar,nH,1);
Hstar = Hmin + (Hmax-Hmin)*rand(nH,1);
Vstar = Vmin + (Vmax - Vmin)*rand(nV,1);
mStar = (betaStarM.*phiStarM.*M)'*Hstar;
rStar = (phiStarM.*M)*Vstar./(1 - a*Hstar/K);
delta = 0.1

r = rStar.*(1 - delta + 2*delta*rand(nH,1));
m =  mStar.*(1 - delta + 2*delta*rand(nH,1));
phi = phiStar.*(1 - delta + 2*delta*rand(1,nV)); 
beta = betaStar.*(1 - delta + 2*delta*rand(1,nV));

phi = repmat(phi,nH,1);
beta = repmat(beta,nH,1);

xeq = equilibrium(M,r,a,K,phi,beta,m);
x0  = abs((0.95 +  2*rand(size(xeq))).*xeq );

% if rank(M) == 10
%     xeq = equilibrium(M,r,a,K,phi,beta,m)
%     x0  = (0.95 +  2*rand(size(xeq))).*xeq ;
% else
%     x0 = rand(nH+nV,1)
% end

tic
[t,x] = predator_prey_integrator(M,r,a,K,phi,beta,m,x0,t);  %x = [H;V]
toc

% Plots
figure
setfigure(12,10,60,10)
semilogy(t,x)
figure
setfigure(12,10,72,10)
plot(t,x)

name = ['data/series_' int2str(nH) '_' int2str(nV) '_2']
save(name)


