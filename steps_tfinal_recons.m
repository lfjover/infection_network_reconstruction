clear all
load('data/params/params_feasible')
nExp = 1;
tfinalV = 10:4:100;
stepsV = 10:10:200;
matV = 1:100;
[nH,nV,nM] = size(matrices);

dt = 0.01;

tVector = dt:dt:tfinalV(end);

reconsErrorA = zeros(length(tfinalV), length(stepsV), length(matV));
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

    runName = ['iM_' num2str(iM)];
    timeSeries = ['data/tseries/tseries_steps_tfinal'];
    save(timeSeries,'allRuns','dt','M')

    
    i_tfinal = 0;
    for tfinal = tfinalV
        i_tfinal = i_tfinal + 1;
        i_steps = 0;
        for steps = stepsV
            i_steps = i_steps + 1;
            nMeas = floor(tfinal/dt/steps);
            lambda1 = 0;
            lambda2 = 0;
            [W,H,Mrec,mrec] = fun_net_recons(timeSeries,[1:nExp],...
                                             nMeas,steps,lambda1,lambda2);
            reconsErrorA(i_tfinal,i_steps,i_mat) = error_recons(Mrec,Mtilde);
            ['matrix: ' num2str(i_mat) ', tfinal:' num2str(i_tfinal) ...
                ', steps:' num2str(i_steps)]
        end
    end
end
save('data/rec_steps_tfinal_2','reconsErrorA','tfinalV','stepsV','matV','dt')
runTime = toc
%% figure
close all
load('data/rec_steps_tfinal.mat')
reconsErrorM = mean(reconsErrorA,3);

figure
width = 25;
height = 17;
fs = 20;
setfigure(width,height,70,16)

% error
imagesc(flipud(reconsErrorM))
hold on
colormap jet
colorbar

%isoclines
N_vector = [100 200 400];
dt_iso = [0; 0.21];%rescale for the new x axis determined with set
x_iso_mat = repmat(dt_iso*100,1,length(N_vector));
y_iso_mat = [];
for N = N_vector;
    'heyt'
    T = dt_iso*N;
    y_iso = 23 - (T-10)./4; %rescale for the new x axis determined with set
    y_iso_mat = [y_iso_mat ,y_iso];
end
plot(x_iso_mat, y_iso_mat, '-k', 'linewidth', 3)
hold off

% % text for isoclines 100 200 400
% x = 0.08;
% x_text = x*100;
% N = 400;
% y_text = 23 - (x*N-10)./4;
% 
% text(x_text+0.5, y_text, {'1000',  'meas.'},'fontweight', 'bold', ...
%         'color', 'white','fontsize', 15)
%     
% x = 0.13;
% x_text = x*100;
% N = 500;
% y_text = 23 - (x*N-10)./4;
% 
% text(x_text, y_text+1, {'500',  'meas.'},'fontweight', 'bold', ...
%         'color', 'black','fontsize', 15)
%     
% x = 0.16;
% x_text = x*100;
% N = 250;
% y_text = 23 - (x*N-10)./4;
% 
% text(x_text, y_text+1, {'250',  'meas.'},'fontweight', 'bold', ...
%         'color', 'black','fontsize', 15)

% text for isoclines 250 500 1000
x = 0.17;
x_text = x*100;
N = 400;
y_text = 23 - (x*N-10)./4;

text(x_text, y_text-2.5, {'400',  'meas.'},'fontweight', 'bold', ...
        'color', 'black','fontsize', 15)
    
x = 0.17;
x_text = x*100;
N = 200;
y_text = 23 - (x*N-10)./4;

text(x_text, y_text-1.7, {'200',  'meas.'},'fontweight', 'bold', ...
        'color', 'black','fontsize', 15)
    
x = 0.17;
x_text = x*100;
N = 100;
y_text = 23 - (x*N-10)./4;

text(x_text, y_text-1.3, {'100',  'meas.'},'fontweight', 'bold', ...
        'color', 'black','fontsize', 15)


set(gca,'ytick',1:length(tfinalV),'yticklabel',fliplr(tfinalV),...
       'xtick',1:length(stepsV),'xticklabel',stepsV*dt)
xlabel('$\Delta t$ (hours)', 'interpreter', 'latex', 'fontsize', fs)
ylabel('Total hours', 'interpreter', 'latex', 'fontsize', fs)
title('$\overline{Error}_{rec}$',...
        'interpreter','latex','fontsize',fs+5)
print('-dpdf','../manuscript/figures/error_deltat_finalt.pdf')