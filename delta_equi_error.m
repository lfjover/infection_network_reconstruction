% average error vs delta from equilibrium, an ensemble of 100 matrices
clear all
load('data/params/params_feasible')

steps = 10;
tfinal = 96;
matV = 1:20:100;
deltaEquiV = 0.01:0.01:0.5;
[nH,nV,nM] = size(matrices);

dt = 0.01;
tVector = dt:dt:tfinal;

reconsErrorM = zeros(length(matV), length(deltaEquiV));
%condNumM = zeros(length(matV), length(deltaEquiV));
%cvHM = zeros(length(matV), length(deltaEquiV));
i_mat = 0;
for iM = matV
    iSet = iM;
    i_mat = i_mat + 1;
    M = matrices(:,:,iM);
    [r,a,phi,beta,m] = params{iSet,:}; 
    Mtilde = M.*phi.*beta;
    
    i_delta = 0;
    for delta = deltaEquiV
         i_delta = i_delta + 1;
        xEq = equilibrium(M,r,a,K,phi,beta,m);
        x0 = xEq.*((1-delta) + (2*delta)*rand(nV+nH,1));

        allRuns = cell(1,2);
        [t,x] = predator_prey_integrator(M,r,a,K,phi,...
                beta,m,x0,tVector); 
        allRuns(1,:) = {t,x};


%         runName = ['iM_' num2str(iM) '_delta_' num2str(i_delta)];
%         timeSeries = ['data/tseries/tseries_single_' runName];
%         save(timeSeries,'allRuns','dt','M')
        timeSeries = 'timeSeries';
        save(timeSeries,'allRuns','dt','M')
    
        nMeas = floor(tfinal/dt/steps);
        lambda1 = 0;
        lambda2 = 0;
        [W,H,Mrec,mrec] = fun_net_recons(timeSeries,[1],...
                                         nMeas,steps,lambda1,lambda2);
        reconsErrorM(i_mat,i_delta) = error_recons(Mrec,Mtilde);
        %svdH = svd(H);
        %condNumM(i_mat,i_delta) = max(svdH)/min(svdH); 
        %cvHM(i_mat,i_delta) = mean(std(x(:,1:nH))./mean(x(:,1:nH)));

    end
end
%save('data/rec_delta_equi','reconsErrorM','matV','deltaEquiV')

%% mean error vs delta from equilibrium
%load('data/rec_delta_equi.mat')

figure
fs = 20;
plot(deltaEquiV, mean(reconsErrorM), '-ob', 'linewidth',3)

legend boxoff
xlabel('$\delta$ from equilibrium','fontsize',fs,'interpreter','latex')
ylabel('$\overline{Error}_{rec}$',...
        'interpreter','latex','fontsize',fs)
setfigure(15,10,68,6)
% hold on
% error_std = std(reconsErrorM);
% plot(deltaEquiV, mean(reconsErrorM)+error_std, '--k', 'linewidth',3)
% plot(deltaEquiV, mean(reconsErrorM)-error_std, '--k', 'linewidth',3)