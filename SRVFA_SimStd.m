function [simStd,bestAngles,sigparams] = SRVFA_SimStd(M0,TR,T1,T2star,FA1,FA2,T1dovT1,zeta,TE,sigma,aveDims,doT2starCor,doSRVFA_T2starAll,varargin)
% SRVFA_SimStd: simulates the standard deviation of the T1 measurement
% using the single-reference variable flip angle method (Svedin et al. 2019)
% from the analytic variance equation presented by Malmberg et al. 2021
% based on the spoiled gradient echo signal parameters and other parameters
% described by Malmberg et al. The T2* correction can also be simulated by
% using multiple TEs at a time in this call.

% Inputs
if ~exist('doT2starCor','var') || isempty(doT2starCor)
    doT2starCor = 0;
end
if ~exist('doSRVFA_T2starAll','var') || isempty(doT2starCor) || doT2starCor == 0
    doSRVFA_T2starAll = 0;
end

if ~isempty(varargin)
    baselineAvg = varargin{1};
    dynamicAvg = varargin{2};
    try 
        if mod(baselineAvg,1) ~=0 || mod(dynamicAvg,1) ~=0 || baselineAvg < 1 || dynamicAvg < 1
            error('baseline and dynamic averages must be positive integers');
        end
    catch
        error('baseline and dynamic averages must be positive integers');
    end
else
    baselineAvg = 1; dynamicAvg = 1;
end
if ~exist('aveDims','var')
    aveDims = [];
end
valids = ["TR","T1","T2star","T1dovT1","zeta","TE"];
validsNum = [2 3 4 7 8 9]; 
for ii = 1:length(aveDims)
    if isnumeric(aveDims(ii)) || ~ismember(aveDims(ii),valids)
        disp(valids);
        error('aveDims can only include the above listed items');
    end
end
if ~isempty(aveDims)
    aveDims = validsNum((ismember(valids,aveDims)));
end
tmp = [length(M0) length(TR) length(T1) length(T2star) length(FA1) length(FA2) length(T1dovT1) length(zeta) length(TE)];
tmp2 = tmp(1:end-1); % don't include TE for this one
FA1t = FA1; FA2t = FA2;
[sigs,sigparams] = makeSPGRSigs(M0,TR,T1,TE,T2star,FA1,FA2,T1dovT1,zeta,sigma);
if size(sigparams,2) > 9
    sigma = sigparams(:,10);
else
    sigma = sigma(1); % just use the first one if the TR and sigma dimensions don't match up.
end
nTE = length(TE);
FA1 = sigparams(:,5); FA1 = reshape(FA1,[],nTE);
FA2 = sigparams(:,6); FA2 = reshape(FA2,[],nTE);
FAs = [FA1(:,1),FA2(:,1)]; % only need the first echo's FAs. They are the same for all
TEv = reshape(TE,1,nTE);
TRs = sigparams(:,2); TRs = reshape(TRs,[],nTE); 
TRs = TRs(:,1); % only keep the first echo's TR. They are the same for all.
T1 = reshape(T1,1,1,[]); T1dovT1 = reshape(T1dovT1,1,1,1,1,1,1,[]);
T1ds = T1.*T1dovT1;

sigs1 = sigs(:,1); sigs2  = sigs(:,2); sigs3 = sigs(:,3);
% Now make the TE dimension the second dimension, which is the last one in
% this case

sigs1 = reshape(sigs1,[],nTE); 
sigs2 = reshape(sigs2,[],nTE); 
sigs3 = reshape(sigs3,[],nTE);
% Now I can feed these into the variance calculator I made

if doSRVFA_T2starAll == 1
    variance = SRVFAvarT1dAll(sigs1,sigs2,sigs3,FAs,TEv,TRs,sigma,doT2starCor,baselineAvg,dynamicAvg);
    simVar = (reshape(variance,tmp2));
else
    variance = SRVFAvarT1d(sigs1,sigs2,sigs3,FAs,TEv,TRs,sigma,doT2starCor,baselineAvg,dynamicAvg);
    simVar = (reshape(variance,tmp));
end

if ~isempty(aveDims) 
    if isnumeric(aveDims)
        simStd = mean(sqrt(simVar)./T1ds,aveDims); % the stdev should be compared relative to T1d, not the variance
    end
else
    simStd = sqrt(simVar)./T1ds;
end

if nargout > 1
    [grid1,grid2] = ndgrid(FA1t,FA2t);
    mask = false(1,1,1,1,length(FA1t),length(FA2t));
    mask(grid1>= (0.9.*grid2)) = true;
    tmp2 = simStd;
    tmp2(mask) = Inf;
    [~,bestInds] = min(tmp2,[],[5 6],'linear');
    [~,~,~,~,bestRefInd,bestDynInd,~,~,~,~,~] = ind2sub(size(tmp2),bestInds);
%     [~,bestDynInd] = min(tmp3,[],6);
%     [~,bestRefInd] = min(tmp2(:,:,:,:,:,bestDynInd,:,:,:,:,:),[],5);
    lastDim = length(size(bestRefInd));
    bestAngles = cat(lastDim+1,FA1t(bestRefInd), FA2t(bestDynInd));
end


