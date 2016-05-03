function [W,H,MpbRec,mrec] = fun_net_recons(timeSeries,...
                                expV,nMeas,steps,lambda1,lambda2)
% reconstruct host and virus densities
% from multiple experiments, each one with all of the host present.
%
%timeSeries: .mat file with time series
%expV: vector with experiments to include
%nMeas: measurements per experiment ( rows in Hm)
%noiseSd: standard deviation of the noise
%steps: diff between consecutive time measurements: 1 time step is 0.01 hr
%lambda: L1 regularization coeficient 

d = load(timeSeries);
%Mtilde = d.M.*d.phi.*d.beta;
W = [];
H = [];
[nH,nV] = size(d.M);

for expe = expV
    [~,x] = d.allRuns{expe,:}; %data from experiment expe
    v = x(:,nH+1:end); %virus densities,complete time series, each colunm is one virus
    h = x(:,1:nH); % host densities
    Wexp = [];
    Hexp = [];
    % Sampling the time series 

    for inIdx = 1:steps:nMeas*steps
        vInit = v(inIdx,:)';
        vFinal = v(inIdx+steps,:)';
        Wexp = [Wexp, log(vFinal./vInit)/(d.dt*steps)];
        Hexp = [Hexp, h(inIdx,:)'];
    end
    W = [W, Wexp];
    H = [H, Hexp];
end

% reconstruct infection matrix
cvx_begin
    variables MpbRecT(nV,nH) mrec(nV,1)
    minimize( norm([MpbRecT,-mrec]*[H;ones(1,nMeas*length(expV))]-W,'fro')...
        + lambda1*norm(MpbRecT(:),2) ...
        + lambda2*norm(mrec,2))
    subject to
    0 <= MpbRecT
    0 <= mrec
cvx_end

MpbRec = MpbRecT.';

