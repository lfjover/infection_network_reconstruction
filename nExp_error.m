clear all
load('data/params/params_feasible')

nExpMax= 40;
delta = 0.2;
steps = 10;
tfinal = 96;
matV = 1:100;
deltaEquiV = 0.01:0.01:0.5;
[nH,nV,nM] = size(matrices);

dt = 0.01;
tVector = dt:dt:tfinal;

reconsErrorM = zeros(length(matV), length(deltaEquiV));
condNumM = zeros(length(matV), length(deltaEquiV));
cvHM = zeros(length(matV), length(deltaEquiV));
i_mat = 0;
for iM = matV
    iSet = iM;
    i_mat = i_mat + 1;
    M = matrices(:,:,iM);
    [r,a,phi,beta,m] = params{iSet,:}; 
    Mtilde = M.*phi.*beta;
    
    i_delta = 0;
    for expe = 1:nExpMax
        i_delta = i_delta + 1;
        xEq = equilibrium(M,r,a,K,phi,beta,m);
        x0 = xEq.*((1-delta) + (2*delta)*rand(nV+nH,1));

        allRuns = cell(nExp,2);
        [t,x] = predator_prey_integrator(M,r,a,K,phi,...
                beta,m,x0,tVector); 
        allRuns(1,:) = {t,x};
    end

        runName = ['iM_' num2str(iM) '_delta_' num2str(i_delta)];
        timeSeries = ['data/tseries/tseries_single_' runName];
        save(timeSeries,'allRuns','dt','M')

    for nExp = 1:nExpMax
        nMeas = floor(tfinal/dt/steps/nExp);
        lambda1 = 0;
        lambda2 = 0;
        [W,H,Mrec,mrec] = fun_net_recons(timeSeries,[1:nExp],...
                                         nMeas,steps,lambda1,lambda2);
        reconsErrorM(i_mat,i_delta) = error_recons(Mrec,Mtilde);
        svdH = svd(H);
        condNumM(i_mat,i_delta) = max(svdH)/min(svdH); 
        cvHM(i_mat,i_delta) = mean(std(x(:,1:nH))./mean(x(:,1:nH)));
    end
end
save('data/rec_nExp','reconsErrorM','matV','deltaEquiV')

%%
load('data/rec_nExp.mat')

figure
width = 10;
height = 8;
fs = 12; 
setfigure(width,height,70,16)

plot(deltaEquiV,mean(reconsErrorM),'o-')

xlabel('$\delta$', 'interpreter', 'latex', 'fontsize', fs)
ylabel('$Error_{rec}$', 'interpreter', 'latex', 'fontsize', fs)
%%
figure
width = 10;
height = 8;
fs = 12;
setfigure(width,height,70,16)

semilogx(condNumM(:),reconsErrorM(:),'o')

xlabel('condNum', 'interpreter', 'latex', 'fontsize', fs)
ylabel('$Error_{rec}$', 'interpreter', 'latex', 'fontsize', fs)
%%
figure
width = 10;
height = 8;
fs = 12;
setfigure(width,height,70,16)

plot(cvHM(:),reconsErrorM(:),'o')
ylim([0 1])

xlabel('cvHM', 'interpreter', 'latex', 'fontsize', fs)
ylabel('$Error_{rec}$', 'interpreter', 'latex', 'fontsize', fs)