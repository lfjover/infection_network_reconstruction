clear all
load('data/params/params_feasible')

nExp = 20;
steps = 10;
tfinal = 96;
matV = 1:100;
snrV = 1:40;
dt = 0.01;

nMeasExp = floor(tfinal/dt/steps/nExp);
lambda1 = 0;
lambda2 = 0;

reconsErrorM = zeros(length(matV),length(snrV));

i_mat = 0;
tic
for iM = matV
    iSet = iM;
    i_mat = i_mat + 1;
    M = matrices(:,:,iM);
    [r,a,phi,beta,m] = params{iSet,:}; 
    Mtilde = M.*phi.*beta;

    runName = ['iM_' num2str(iM)];
    timeSeries = ['data/tseries/tseries_multi_' runName];
    load(timeSeries)
    
    i_snr = 0;
    for snrDB = snrV      
        i_snr = i_snr +1 ;
        [ 'matrix: ' num2str(iM) ', i_snr: ' num2str(i_snr)]
            
        [W,H,Mrec,mrec] = fun_net_recons_snr(timeSeries,[1:nExp],...
                                         nMeasExp,steps,lambda1,lambda2,snrDB);
        reconsErrorM(i_mat,i_snr) = error_recons(Mrec,Mtilde);
    end
    
    
end
runTime = toc
save('data/rec_snr','reconsErrorM')

%% Reconstruction error vs noise
load('data/rec_snr')
fs = 15;
figure
plot(snrV, mean(reconsErrorM),'-b', 'markersize',8,'linewidth',3)
hold on
plot(snrV, mean(reconsErrorM)+std(reconsErrorM),'--k', 'markersize',8,'linewidth',3)
plot(snrV, mean(reconsErrorM)-std(reconsErrorM),'--k', 'markersize',8,'linewidth',3)
hold off
xlabel('SNR (dB)','fontsize',fs)
ylabel('$\overline{Error}_{rec}$','interpreter','latex','fontsize',fs)
setfigure(10,7,68,6)
