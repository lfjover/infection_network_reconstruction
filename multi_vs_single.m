clear all
load('data/params/params_feasible')

nExpMax= 20;
nExpV = 1:nExpMax;
steps = 10;
tfinal = 96;
matV = 1:100;
delta= 0.2;
[nH,nV,nM] = size(matrices);

dt = 0.01;
tVector = dt:dt:tfinal;

reconsErrorSingle = zeros(length(nExpV),length(matV));
reconsErrorMulti = zeros(1,length(matV));

allRuns = cell(nExpMax,2);
i_mat = 0;
tic
for iM = matV
    iSet = iM;
    i_mat = i_mat + 1;
    M = matrices(:,:,iM);
    [r,a,phi,beta,m] = params{iSet,:}; 
    Mtilde = M.*phi.*beta;
    
    allRuns = cell(nExpMax,2);
    for expe = nExpV
        xEq = equilibrium(M,r,a,K,phi,beta,m);
        x0 = xEq.*((1-delta) + (2*delta)*rand(nV+nH,1));

        [t,x] = predator_prey_integrator(M,r,a,K,phi,...
                beta,m,x0,tVector); 
        allRuns(expe,:) = {t,x};
    end

    runName = ['iM_' num2str(iM)];
    timeSeries = ['data/tseries/tseries_multi_' runName];
    save(timeSeries,'allRuns','dt','M')

    %reconstruction multiple experiment 
    nMeas = floor(tfinal/dt/steps/nExpMax);
    lambda1 = 0;
    lambda2 = 0;
    [W,H,Mrec,mrec] = fun_net_recons(timeSeries,nExpV,...
                                     nMeas,steps,lambda1,lambda2);
    reconsErrorMulti(i_mat) = error_recons(Mrec,Mtilde);
    figure;
    plot_recons(Mtilde,Mrec,H,W,0)
    drawnow
    
    i_expe = 0;
    for expe = nExpV
        i_expe = i_expe +1 ;
        nMeasExp = floor(tfinal/dt/steps);
        lambda1 = 0;
        lambda2 = 0;
        [W,H,Mrec,mrec] = fun_net_recons(timeSeries,expe,...
                                         nMeasExp,steps,lambda1,lambda2);
        if expe < 5
         figure;
        plot_recons(Mtilde,Mrec,H,W,0)
        drawnow
        end
        reconsErrorSingle(i_expe,i_mat) = error_recons(Mrec,Mtilde);
    end
    
    
end
runTime = toc
save('data/rec_multi_v_single','reconsErrorMulti','reconsErrorSingle',...
       'matV')
%%
%load('data/rec_multi_v_single.mat')
 load('data/params/matrices_10by10_100_invertible.mat')
figure;
setfigure(25,10,68,6)
fs = 20;
plot(sortNest, reconsErrorMulti, '-ok', 'linewidth',3)
hold on
plot(sortNest, mean(reconsErrorSingle), '-ob', 'linewidth',3)
hold off
% hold on
% plot(matV, mean(reconsErrorSingle)+std(reconsErrorSingle),...
%     '--k', 'markersize',8,'linewidth',3)
% plot(matV, mean(reconsErrorSingle)-std(reconsErrorSingle),...
%     '--k', 'markersize',8,'linewidth',3)
% hold off
legend boxoff
xlabel('Nestedness','fontsize',fs,'interpreter','latex')
ylabel('$\overline{Error}_{rec}$',...
        'interpreter','latex','fontsize',fs)
