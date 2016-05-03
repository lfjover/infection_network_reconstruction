% good reconstruction
clear all
close all
load('data/params/first_example')  % load params with model parameters and matrices with matrix ensemble
[r,a,phi,beta,m] = params{1,:};  % model parameters 
[nH,nV] = size(M); % infection matrix is M
Mtilde = M.*phi.*beta;

tfinal = 240;
dt = 0.01;
tVector = dt:dt:tfinal;

% time series are saved in the cell allRuns ( this case is only one time series but is easier to
% use the cell struct for generalizing the code to multiple experiments with the same parameters)
nExp = 1; 
allRuns = cell(nExp,2);
 
xEq = equilibrium(M,r,a,K,phi,beta,m); % calculates equilibrium densities
delta = 0.1; 
x0 = xEq.*((1-delta) + (2*delta)*rand(nV+nH,1));

%run dynamics
[t,x] = predator_prey_integrator(M,r,a,K,phi,...
        beta,m,x0,tVector); 

allRuns(1,:) = {t,x};

%timeSeries = 'data/tseries/tseries_first_example';
timeSeries = 'data/tseries/tseries_example';
save(timeSeries,'allRuns','dt','M')

% calculate reconstruction
steps = 10; %number of steps dt between measurements
nMeas = floor(tfinal/dt/steps);
lambda1 = 0; % not used in the analysis (regularization factors)
lambda2 = 0;
[W,H,Mrec,mrec] = fun_net_recons(timeSeries,1,...
                                             nMeas,steps,lambda1,lambda2);
reconsError = error_recons(Mrec,Mtilde)
%plot_recons(Mtilde,Mrec,H,W,0)
%save('data/rec_first_example','W','H','Mrec','mrec','steps','nMeas','Mtilde')

% figure
plot_recons(Mtilde, Mrec, H, W, timeSeries)