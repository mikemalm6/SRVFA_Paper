function [variance] = SRVFAvarT1d(sigs1,sigs2,sigs3, FAs,TEv,TR,sigma,doT2starCor,varargin)

if ~isempty(varargin)
    if length(varargin) >= 1
        baselineAvg = varargin{1}; 
    end
    if length(varargin) == 2
        dynamicAvg = varargin{2};
    end
else
    baselineAvg = 1;
    dynamicAvg = 1;
end
if ~exist('doT2starCor','var') || isempty(doT2starCor)
    doT2starCor = 0;
end
nTE = length(TEv);
if size(sigs1,1) == nTE; sigs1 = sigs1'; end % make the TE dimension the 2nd dimension
if size(sigs2,1) == nTE; sigs2 = sigs2'; end % make the TE dimension the 2nd dimension
if size(sigs3,1) == nTE; sigs3 = sigs3'; end % make the TE dimension the 2nd dimension
if size(TEv,1) == nTE; TEv = TEv'; end % make the TE dimension the 2nd dimension
if any(size(sigma) == 1)
    if length(sigma) ~= size(sigs1,1) && ~isscalar(sigma)
        sigma = reshape(sigma,[],nTE); sigma = sigma(:,1); % sigma for all echoes will be the same
    end
end

[nR1, nC1] = size(sigs1); [nR2, nC2] = size(sigs2); [nR3, nC3] = size(sigs3); 

if nTE == 1 || (nC1 == 1 && nC2 == 1 && nC3 == 1) || doT2starCor == 0
    doSRVFAT2 = 0;
else
    if nTE ~=1 && (nC1 == nTE && nC2 == nTE && nC3 == nTE) && doT2starCor ~= 0
        doSRVFAT2 = 1;
    else
        error('number of signals is incompatible with TEs provided');
    end
end
if nR1 ~= nR2 || nR1~= nR3
    error('number of signals must be the same for each type');
else
    nR = nR1;
end
% doSRVFAT2 = 1;
% nR = size(sigs1,1);
rho = sigs2./sigs1;
rho2 = sigs3./sigs1;

if numel(FAs) == 2
    a = FAs(1); 
    b = FAs(2);
else
    a = FAs(:,1); % allow for simulation studies to have multiple flip angles 
    b = FAs(:,2); % allow for simulation studies to have multiple flip angles 
end

sina = sind(a); % (nR x 1 vector) or scalar
sinb = sind(b); % (nR x 1 vector) or scalar
cosa = cosd(a); % (nR x 1 vector) or scalar
cosb = cosd(b); % (nR x 1 vector) or scalar

E1 = (rho.*sina - sinb)./(rho.*sina.*cosb - sinb.*cosa); % (nR x nTE matrix)
E1est = (rho2.*sina - sinb)./(rho2.*sina.*cosb - sinb.*cosa); % (nR x nTE matrix)
clear rho; clear rho2;
gamma = (1-E1)./(1-E1.*cosa).*(1-E1est.*cosa)./(1-E1est.*cosb); % (nR x nTE matrix)

% Now work on dE2rdS
sumTE = sum(TEv);
enddeds1 = (nTE.*TEv./sigs1 - sum(TEv)./sigs1).*(nTE.*sum(TEv.*TEv) - (sumTE.*sumTE))./((nTE.*sum(TEv.*log(sigs1),2) - sum(TEv).*sum(log(sigs1),2)).^2); % (nR x nTE matrix) 
enddeds2 = (nTE.*TEv./sigs2 - sum(TEv)./sigs2).*(nTE.*sum(TEv.*TEv) - (sumTE.*sumTE))./((nTE.*sum(TEv.*log(sigs2),2) - sum(TEv).*sum(log(sigs2),2)).^2); % (nR x nTE matrix)
enddeds3 = (nTE.*TEv./sigs3 - sum(TEv)./sigs3).*(nTE.*sum(TEv.*TEv) - (sumTE.*sumTE))./((nTE.*sum(TEv.*log(sigs3),2) - sum(TEv).*sum(log(sigs3),2)).^2); % (nR x nTE matrix)


if doSRVFAT2 == 1
    [E2rg,T2star,T2stard] = SRVFA_E2r(sigs1,sigs2,sigs3,TEv,TEv); % E2rg = (nR x nTE matrix), T2star = (nR x 1 vector), T2stard = (nR x 1 vector)
else
    E2rg = 1;
end
T1dest = real(-TR./log((E2rg-gamma)./(E2rg-gamma.*cosb))); % (nR x nTE matrix)
lambda = SRVFA_lambda(E2rg,gamma, b, T1dest, TR); % (nR x nTE matrix)
bot21 = (sigs2./sigs1.*sina.*cosb - sinb.*cosa); bot21 = bot21.*bot21;
bot31 = (sigs3./sigs1.*sina.*cosb - sinb.*cosa); bot31 = bot31.*bot31;
df2_1 = (-sigs2./(sigs1.*sigs1)).*sina.*sinb.*(cosb-cosa)./bot21; % (nR x nTE matrix)
df3_1 = (-sigs3./(sigs1.*sigs1)).*sina.*sinb.*(cosb-cosa)./bot31; % (nR x nTE matrix)
df2_2 = (1./sigs1).*sina.*sinb.*(cosb-cosa)./bot21; % (nR x nTE matrix)
df3_3 = (1./sigs1).*sina.*sinb.*(cosb-cosa)./bot31; % (nR x nTE matrix)


dgds = zeros(nR,nTE*3);
dgds(:,1:3:end) = df2_1.*(cosa - 1)./((1-E1.*cosa).^2).*(1-E1est.*cosa)./(1-E1est.*cosb) + ...
   df3_1.*(cosb-cosa)./((1-E1est.*cosb).^2).*(1-E1)./(1-E1.*cosa); % (nR x nTE matrix)
dgds(:,2:3:end) = df2_2.*(cosa - 1)./((1-E1.*cosa).^2).*(1-E1est.*cosa)./(1-E1est.*cosb); % (nR x nTE matrix)
dgds(:,3:3:end) = df3_3.*(cosb-cosa)./((1-E1est.*cosb).^2).*(1-E1)./(1-E1.*cosa); % (nR x nTE matrix)
% dgds is a nR x 3*nTE matrix)


if doSRVFAT2 == 1
    TEv2 = reshape(TEv,1,1,[]);
    beg = TEv2.*exp(-TEv2.*(1./T2stard - 1./T2star)); % (nR x 1 x nTE matrix)
    deds1 = -beg./(2.*T2star.*T2star).*enddeds1; 
    deds2 = -beg./(2.*T2star.*T2star).*enddeds2; 
    deds3 = beg./(T2stard.*T2stard).*enddeds3; 
    deds = zeros(nR,3*nTE,nTE);
    deds(:,1:3:end,:) = deds1; clear deds1; 
    deds(:,2:3:end,:) = deds2; clear deds2;
    deds(:,3:3:end,:) = deds3; clear deds3;
end
avgs = zeros(1,3*nTE);
avgs(1:3:end) = baselineAvg;
avgs(2:3:end) = baselineAvg;
avgs(3:3:end) = dynamicAvg;
variance = zeros(nR,nTE);
for ii = 1:nTE
    userange = (1:3) + (3*(ii-1));
    mask = zeros(1,3*nTE); mask(userange) = 1;
%     mask = ismember((1:(3*nTE)),userange);
    mask = repmat(mask,[nR,1]);
    dgdsuse = dgds.*mask; % (nR x 3nTE matrix)
    if doSRVFAT2 == 1
        dedsuse = deds(:,:,ii); % (nR x 3NTE vector)
        tmpvarend = (E2rg(:,ii).*dgdsuse - gamma(:,ii).*dedsuse); tmpvarend = tmpvarend .* tmpvarend;
        varend = sum(1./avgs.*tmpvarend,2);
        variance(:,ii) = sigma.*sigma.*(lambda(:,ii).*lambda(:,ii)).*varend; % should be vector of length TEv
    else
        tmpvarend = (E2rg.*dgdsuse); tmpvarend = tmpvarend.*tmpvarend;
        varend = sum(1./avgs.*tmpvarend,2);
        variance(:,ii) = sigma.*sigma.*(lambda(:,ii).*lambda(:,ii)).*varend; % should be vector of length TEv
    end
end 

end
    

function lambda = SRVFA_lambda(E2rg,gamma, betaAng, T1d_est, TR)
cosb = cosd(betaAng);
lambda = T1d_est.^2.*(cosb-1)./(TR.*(E2rg-gamma).*(E2rg-gamma.*cosb));
end

function [E2r,T2star,T2stard,T2star1,T2star2] = SRVFA_E2r(sigs1,sigs2,sigs3,TEv,TE)
nTE = length(TEv);
TE = reshape(TE,1,[]);
TEv = reshape(TEv,1,[]);
top = -(nTE.*sum(TEv.^2) - (sum(TEv).^2)); % scalar
bot1 = nTE.*sum(TEv.*log(sigs1),2) - sum(TEv).*sum(log(sigs1),2); % (nR x 1 vector)
bot2 = nTE.*sum(TEv.*log(sigs2),2) - sum(TEv).*sum(log(sigs2),2); % (nR x 1 vector)
bot3 = nTE.*sum(TEv.*log(sigs3),2) - sum(TEv).*sum(log(sigs3),2); % (nR x 1 vector)

T2star1 = top./bot1; T2star2 = top./bot2; % (nR x 1 vector)
T2star = 0.5.*T2star1 + 0.5.*T2star2; % (nR x 1 vector)
T2stard = top./bot3; % (nR x 1 vector)

E2r = exp(-TE.*(1./T2stard - 1./T2star)); % (nR x length(TE) matrix)
end
