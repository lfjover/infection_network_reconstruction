clear all
load('data/params/params_feasible')
nExp = 1;
stepsV = ceil(logspace(0,log10(2400),30));
nMeasV = [100 200 400];
matV = 1:100;
[nH,nV,nM] = size(matrices);

dt = 0.01;

tVector = dt:dt:(nMeasV(end)*stepsV(end)*dt);

reconsErrorA = zeros(length(nMeasV), length(stepsV), length(matV));
i_mat = 0;
tic
for iM = matV
    iSet = iM;
    i_mat = i_mat + 1;
    M = matrices(:,:,iM);
    [r,a,phi,beta,m] = params{iSet,:}; 
    Mtilde = M.*phi.*beta;

    xEq = equilibrium(M,r,a,K,phi,beta,m);
    delta = 0.1;
    x0 = xEq.*((1-delta) + (2*delta)*rand(nV+nH,1));

    allRuns = cell(nExp,2);
    [t,x] = predator_prey_integrator(M,r,a,K,phi,...
            beta,m,x0,tVector); 
    allRuns(1,:) = {t,x};

    %runName = ['iM_' num2str(iM)];
    timeSeries = ['data/tseries/tseries_deltat'];
    save(timeSeries,'allRuns','dt','M')
   
    i_nMeas = 0;
    for nMeas = nMeasV
        i_nMeas = i_nMeas + 1;
        i_steps = 0;
        for steps = stepsV
            i_steps = i_steps + 1;
            lambda1 = 0;
            lambda2 = 0;
            [W,H,Mrec,mrec] = fun_net_recons(timeSeries,[1:nExp],...
                                             nMeas,steps,lambda1,lambda2);
            reconsErrorA(i_nMeas,i_steps,i_mat) = error_recons(Mrec,Mtilde);
            ['matrix: ' num2str(i_mat) ', tfinal:' num2str(i_nMeas) ...
                ', steps:' num2str(i_steps)]
        end
    end
end
save('data/rec_deltat','reconsErrorA','nMeasV','stepsV','matV','dt')
runTime = toc
%% figure
close all
load('data/rec_deltat.mat')
reconsErrorM = mean(reconsErrorA,3);

figure
width = 18;
height = 11;
fs = 22;
setfigure(width,height,70,16)
semilogx(stepsV(10:end)*dt,reconsErrorM(:,10:end),'o-','linewidth',3)
legend({'100 Meas.','200 Meas.','400 Meas.'})
legend boxoff
xlabel('$\Delta t$ (hours)','fontsize',fs,'interpreter','latex')
ylabel('$\overline{Error}_{rec}$',...
        'interpreter','latex','fontsize',fs)
print('-dpdf','../manuscript/submission_rsopen/figures/error_deltat_fixed_nMeas.pdf')
%%
figure
width = 18;
height = 11;
fs = 20;
setfigure(width,height,70,16)
plot(t,x,'linewidth',2)
xlabel('time (h)')
ylabel('Density (particle/ml)')
print('-dpdf','example_dynamics.pdf')
