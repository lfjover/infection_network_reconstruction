%% feasible parameter sets for mod-nestedness ensemble
clear all
load('data/params/matrices_10by10_100_invertible.mat')

[nH,nV,nM] = size(matrices);
matV = 1:nM;

a = ones(nH);
phiMin = 10^-8; phiMax = 10^-7;
betaMin = 10; betaMax = 50;
Hmin = 10^3; Hmax = 10^4;
Vmin = 10^6; Vmax = 10^7;
K = Hmax*100;

params = cell(nM,5);

for iM = matV
    M = matrices(:,:,iM);
    phi = phiMin + (phiMax - phiMin)*rand(nH,nV);
    beta = betaMin + (betaMax - betaMin)*rand(nH,nV);
    Hstar = Hmin + (Hmax-Hmin)*rand(nH,1);
    Vstar = Vmin + (Vmax - Vmin)*rand(nV,1);
    m = (beta.*phi.*M)'*Hstar;
    r = (phi.*M)*Vstar./(1 - a*Hstar/K);
    params(iM,:) = {r,a,phi,beta,m}; 
end

save('data/params/params_feasible','params','matrices','K')


%% first example parameters and matrix
clear all
nH = 10;
nV = 10;
nConn = 10;

%create random matrix with diagonal elements
M = eye(nH,nV);  %start with diagonal matrix
idx_off_diag = reshape(1:nH*nV,nH,nV);
idx_off_diag= idx_off_diag - diag(diag(idx_off_diag));
idx_off_diag = idx_off_diag(:);
idx_off_diag(idx_off_diag ==0) = [];
M(randsample(idx_off_diag,nConn)) = 1;

a = ones(nH);
phiMin = 10^-8; phiMax = 10^-7;
betaMin = 10; betaMax = 50;
Hmin = 10^3; Hmax = 10^4;
Vmin = 10^6; Vmax = 10^7;
K = Hmax*100;

params = cell(1,5);


phi = phiMin + (phiMax - phiMin)*rand(nH,nV);
beta = betaMin + (betaMax - betaMin)*rand(nH,nV);
Hstar = Hmin + (Hmax-Hmin)*rand(nH,1);
Vstar = Vmin + (Vmax - Vmin)*rand(nV,1);
m = (beta.*phi.*M)'*Hstar;
r = (phi.*M)*Vstar./(1 - a*Hstar/K);
params(1,:) = {r,a,phi,beta,m}; 


imagesc(M)

%save('data/params/first_example','params','M','K')