function [varianceAll, variance] = SRVFAvarT1dAll(sigs1,sigs2,sigs3, FAs,TEv,TR,sigma,doT2starCor,varargin)

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
global nTE; 
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
    if doT2starCor == 1
        error('Attempt to use SR-VFA-T2* method with only one echo time invalid');
    end
else
    if nTE ~=1 && (nC1 == nTE && nC2 == nTE && nC3 == nTE) && doT2starCor ~= 0
        doSRVFAT2 = 1;
    else
        error('number of signals is incompatible with TEs provided');
    end
end
if nR1 ~= nR2 || nR1~= nR3
    error('number of signals must be the same for each type');
end

global nR; 
nR = nR1;

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
lambda = SRVFA_lambda(E2rg,gamma, cosb, T1dest, TR); % (nR x nTE matrix)
bot21 = (sigs2./sigs1.*sina.*cosb - sinb.*cosa); bot21 = bot21.*bot21;
bot31 = (sigs3./sigs1.*sina.*cosb - sinb.*cosa); bot31 = bot31.*bot31;
df2_1 = (-sigs2./(sigs1.*sigs1)).*sina.*sinb.*(cosb-cosa)./bot21; % (nR x nTE matrix)
df3_1 = (-sigs3./(sigs1.*sigs1)).*sina.*sinb.*(cosb-cosa)./bot31; % (nR x nTE matrix)
df2_2 = (1./sigs1).*sina.*sinb.*(cosb-cosa)./bot21; % (nR x nTE matrix)
df3_3 = (1./sigs1).*sina.*sinb.*(cosb-cosa)./bot31; % (nR x nTE matrix)
clear bot21; clear bot31;


dgds = zeros(nR,nTE*3);
dgds(:,1:3:end) = df2_1.*(cosa - 1)./((1-E1.*cosa).^2).*(1-E1est.*cosa)./(1-E1est.*cosb) + ...
   df3_1.*(cosb-cosa)./((1-E1est.*cosb).^2).*(1-E1)./(1-E1.*cosa); % (nR x nTE matrix)
dgds(:,2:3:end) = df2_2.*(cosa - 1)./((1-E1.*cosa).^2).*(1-E1est.*cosa)./(1-E1est.*cosb); % (nR x nTE matrix)
dgds(:,3:3:end) = df3_3.*(cosb-cosa)./((1-E1est.*cosb).^2).*(1-E1)./(1-E1.*cosa); % (nR x nTE matrix)
% dgds is a nR x 3*nTE matrix)
clear df3_3 df2_2 df3_1 df2_1;


if doSRVFAT2 == 1
    TEv2 = reshape(TEv,1,1,[]);
    beg = TEv2.*exp(-TEv2.*(1./T2stard - 1./T2star)); % (nR x 1 x nTE matrix)
    deds1 = -beg./(2.*T2star.*T2star).*enddeds1; clear enddeds1;
    deds2 = -beg./(2.*T2star.*T2star).*enddeds2; clear enddeds2
    deds3 = beg./(T2stard.*T2stard).*enddeds3; clear enddeds3
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
%     mask = repmat(mask,[nR,1]);
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

clear varend tmpvarend
%% Now that the variance is determined for each echo, I have T1dest_k and effectively dT1destdS, and w_k. I need to delineate those out
w_kbot = sum(1./variance,2);
w_k = reshape((1./variance)./w_kbot,[nR, 1, nTE]); % sum along the TE dimension. % These are the weights for the linear combination (nR x 1 x nTE matrix)
E2rguse = beg./TEv2; % now it is just E2d/E2_est_k (nR x 1 x nTE matrix) clear beg;
dT1dk_dS = zeros(nR,3*nTE,nTE);
for kk = 1:nTE
    userange = (1:3) + 3*(kk-1);
    mask = zeros(1,3*nTE); mask(userange) = 1;
    mask = repmat(mask,[nR,1]);
    dgdsuse = dgds.*mask; % (nR x 3nTE matrix)
    dT1dk_dS(:,:,kk) = lambda(:,kk).*(E2rguse(:,1,kk).*dgdsuse - gamma(:,kk).*deds(:,:,kk)); % should be nR x 3*nTE x nTE matrix, because they you need an nR x nTE matrix for each T1d,est_k (for each echo time)
    % dT1dk_dS size is nR, ij, TE
end

T1dest = reshape(T1dest,[nR,1,nTE]);
% nRxnTE to nR x 3*nTE
maskdgds = [];
for ii = 1:nTE
    maskdgds = blkdiag(maskdgds,ones(1,3,1));
end
maskdgds = reshape(maskdgds',[1, 3*nTE,nTE]);
dgdsAll = dgds.*maskdgds; % nR x 3*nTE x nTE

% define dSigmaSqdS
avgsSigSq = reshape(avgs,[1 1 length(avgs)]); % along the pq dimension
dT1dk_dSpq = permute(reshape(permute(dT1dk_dS,[1 3 2]),[nR,nTE,1,3*nTE]),[1 3 4 2]); % now it is nR, 1, pq (copy of ij, but in pq spot), nTE
%gammaAll = reshape(repmat(reshape(gamma,[nR,1,nTE]),[1,3,1]),[nR,3*nTE]);  %%%%%%%%%%%%%%%%%%%%CHECK ME %%%%%%%%%%%%%%%%%
% The last thing I need is dT1dk_dSpqdSij, which requires dKapdSij, E2rk,
% dgdSpq, gk, dEdSpq, and then dEEdSij
dKapdSij = SRVFA_dKapdSij(cosb, reshape(T1dest,[nR, 1,1, nTE]), reshape(dT1dk_dS,[nR,3*nTE,1,nTE]),reshape(E2rguse,[nR,1,1,nTE]),reshape(gamma,[nR,1,1,nTE]),reshape(deds,[nR,3*nTE,1,nTE]),reshape(dgdsAll,[nR,3*nTE,1,nTE]),TR);

dT2stard_dSij = SRVFA_dT2stard_dSij(TEv,sigs3); % nR x 3*nTE
dT2star_dSij = SRVFA_dT2star_dSij(TEv,sigs1,sigs2);% nR x 3*nTE

%%% dEEdSij stuff
% [nR,~,nTE] = size(dgdsAll);
dgdsij = reshape(dgdsAll,[nR,3*nTE,1,nTE]); % nR x 3*nTE x 1 x nTE
dgdspq = reshape(dgdsAll,[nR,1,3*nTE,nTE]); % nR x 1 x 3*nTE x nTE
dedsij = reshape(deds,[nR,3*nTE,1,nTE]);
dedspq = reshape(deds,[nR,1,3*nTE,nTE]);
dgdSijSpq = SRVFA_dgdSijSpq(cosa,cosb,sina,sinb,E1,sigs1,sigs2,sigs3,reshape(E1est,[nR,1,1,nTE]));
dedSijSpq = SRVFA_dedSijSpq(TEv,T2star,T2stard,dT2stard_dSij,dT2star_dSij,sigs1,sigs2,sigs3); 
clearvars T2star T2stard dT2stard_dSij dT2star_dSij sigs1 sigs2 sigs3;
dEEdSij = dedsij.*dgdspq + reshape(E2rguse,[nR,1,1,nTE]).*dgdSijSpq - dgdsij.*dedspq - reshape(gamma,[nR,1,1,nTE]).*dedSijSpq;
%%%
clearvars dgdSijSpq dedSijSpq;

% dEEdSij = SRVFA_dEEdSij(reshape(E2rguse,[nR,1,1,nTE]),dgdsAll,reshape(gamma,[nR,1,1,nTE]),deds,cosa,cosb,sina,sinb,reshape(E1,[nR,1,1,nTE]),reshape(sigs1,[nR,1,1,nTE]),reshape(sigs2,[nR,1,1,nTE]),reshape(sigs3,[nR,1,1,nTE]),reshape(E1est,[nR,1,1,nTE]),reshape(TEv,[1,1,1,nTE]),T2stard(:),T2star(:),dT2stard_dSij,dT2star_dSij);

dT1dk_dSpqdSij = dKapdSij.*(reshape(E2rguse,[nR,1,1,nTE]).*reshape(dgdsAll,[nR,1,3*nTE,nTE]) - reshape(gamma,[nR,1,1,nTE]).*reshape(deds,[nR,1,3*nTE,nTE])) + ...
    reshape(lambda,[nR,1,1,nTE]).*dEEdSij; %%%%%%%%%%%%%%%%%%%%%%%%%
clearvars dgdsAll;
dSigmaSqdS = 2.*sigma.*sigma.*sum(1./avgsSigSq.*dT1dk_dSpq.*dT1dk_dSpqdSij,3); %

dw_dS = -dSigmaSqdS./(reshape(variance.*variance,[nR,1,1,nTE]).*w_kbot) + (sum(dSigmaSqdS./(reshape(variance.*variance,[nR,1,1,nTE])),4))./(reshape(variance,[nR,1,1,nTE]).*w_kbot.*w_kbot);

dT1d_dS = sum(reshape(dw_dS,[nR,3*nTE,nTE]).*T1dest + dT1dk_dS.*w_k,3); % nR x 3*nTE x nTE matrix summed to nR x 3*nTE matrix
varianceAll = sigma.*sigma.*sum(1./avgs.*(dT1d_dS.*dT1d_dS),2); % nR x 1 matrix

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS  %%%%% %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lambda = SRVFA_lambda(E2rg,gamma, cosb, T1d_est, TR)
lambda = T1d_est.*T1d_est.*(cosb-1)./(TR.*(E2rg-gamma).*(E2rg-gamma.*cosb));
end

function [E2r,T2star,T2stard,T2star1,T2star2] = SRVFA_E2r(sigs1,sigs2,sigs3,TEv,TE)
% nTE = length(TEv);
TE = reshape(TE,1,[]);
TEv = reshape(TEv,1,[]);
top = -(nTE.*sum(TEv.*TEv) - (sum(TEv).^2)); % scalar
bot1 = nTE.*sum(TEv.*log(sigs1),2) - sum(TEv).*sum(log(sigs1),2); % (nR x 1 vector)
bot2 = nTE.*sum(TEv.*log(sigs2),2) - sum(TEv).*sum(log(sigs2),2); % (nR x 1 vector)
bot3 = nTE.*sum(TEv.*log(sigs3),2) - sum(TEv).*sum(log(sigs3),2); % (nR x 1 vector)

T2star1 = top./bot1; T2star2 = top./bot2; % (nR x 1 vector)
T2star = 0.5.*T2star1 + 0.5.*T2star2; % (nR x 1 vector)
T2stard = top./bot3; % (nR x 1 vector)

E2r = exp(-TE.*(1./T2stard - 1./T2star)); % (nR x length(TE) matrix)
end

function dKapdSij = SRVFA_dKapdSij(cosb, T1destk, dT1dk_dSij,E2rk,gammak,dE2rdSij,dgammadSij,TR)
E2rmG = E2rk - gammak;
E2rmGcB = E2rk - gammak.*cosb;
dKapdSij = T1destk.*(cosb-1).*(2.*dT1dk_dSij.*(E2rmG).*E2rmGcB - T1destk.*((dE2rdSij-dgammadSij).*E2rmGcB + (dE2rdSij-dgammadSij.*cosb).*E2rmG));
dKapdSij = dKapdSij./(E2rmG.*E2rmG.*TR.*E2rmGcB.*E2rmGcB);

end

function dEEdSij = SRVFA_dEEdSij(E2rk,dgdsAll,gammak,deds,cosa,cosb,sina,sinb,E1,sigs1,sigs2,sigs3,E1estk,TEv,T2stard,T2star,dT2stard_dSij,dT2star_dSij) % still need dedSijSpq and dgdSijSpq
% [nR,~,nTE] = size(dgdsAll);
dgdsij = reshape(dgdsAll,[nR,3*nTE,1,nTE]); % nR x 3*nTE x 1 x nTE
dgdspq = reshape(dgdsAll,[nR,1,3*nTE,nTE]); % nR x 1 x 3*nTE x nTE
dedsij = reshape(deds,[nR,3*nTE,1,nTE]);
dedspq = reshape(deds,[nR,1,3*nTE,nTE]);
dgdSijSpq = SRVFA_dgdSijSpq(cosa,cosb,sina,sinb,E1,sigs1,sigs2,sigs3,E1estk);
dedSijSpq = SRVFA_dedSijSpq(TEv,T2star,T2stard,dT2stard_dSij,dT2star_dSij,sigs1,sigs2,sigs3);

dEEdSij = dedsij.*dgdspq + E2rk.*dgdSijSpq - dgdsij.*dedspq - gammak.*dedSijSpq;
end

function dgdSijSpq = SRVFA_dgdSijSpq(cosa,cosb,sina,sinb,E1,sigs1,sigs2,sigs3,E1estk)
E1 = reshape(E1,[nR,1,1,nTE]);
Ea = 1-E1.*cosa; Eeb = 1 - E1estk.*cosb; Eea = 1-E1estk.*cosa;
%needs to be nR x 3*nTE x nTE;
df21ij = SRVFA_dfyxij(sina,sinb,cosa,cosb,sigs2,sigs1,2);
df21pq = SRVFA_dfyxpq(sina,sinb,cosa,cosb,sigs2,sigs1,2);
df31ij = SRVFA_dfyxij(sina,sinb,cosa,cosb,sigs3,sigs1,3);
df31pq = SRVFA_dfyxpq(sina,sinb,cosa,cosb,sigs3,sigs1,3);

df21_dSpqij = SRVFA_dfyxpqij(sina,sinb,cosa,cosb,sigs2,sigs1,2);
df31_dSpqij = SRVFA_dfyxpqij(sina,sinb,cosa,cosb,sigs3,sigs1,3);
tmp = (cosa-1).*Ea.*(Ea.*df21_dSpqij + 2.*cosa.*df21ij.*df21pq)./(Ea.*Ea.*Ea.*Ea) .* (Eea./Eeb); % term1
tmp = tmp + df31ij.*(cosb-cosa)./(Eeb.*Eeb) .* (df21pq.*(cosa-1)./(Ea.*Ea)); %term1 + term2
tmp = tmp + (cosb-cosa).*Eeb.*(Eeb.*df31_dSpqij + 2.*cosb.*df31ij.*df31pq)./(Eeb.*Eeb.*Eeb.*Eeb) .* (1-E1)./(Ea); %=term1 + term2 + term3
dgdSijSpq = tmp + df21ij.*(cosa-1)./(Ea.*Ea) .* (df31pq.*(cosb-cosa)./(Eeb.*Eeb)); %=term1 + term2 + term3 + term4
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% These are helper functions for SRVFA_dgdSijSpq %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dfyxpqij = SRVFA_dfyxpqij(sina,sinb,cosa,cosb,sigsy,sigsx,y,dsyxdsij,dsyxdspq)
tmp2 = reshape(sigsy./sigsx.*sina.*cosb - sinb.*cosa,[nR,1,1,nTE]);
tmp = reshape(sina.*sinb.*(cosb-cosa)./(tmp2.*tmp2.*tmp2),[nR,1,1,nTE]);
if nargin == 7
    dfyxpqij = tmp.*(tmp2.*SRVFA_dSyxdSpqdSij(sigsy,sigsx,y) - 2.*sina.*cosb.*SRVFA_dSyxdSpq(sigsy,sigsx,y).*SRVFA_dSyxdSij(sigsy,sigsx,y));
elseif nargin == 9
    dfyxpqij = tmp.*(tmp2.*SRVFA_dSyxdSpqdSij(sigsy,sigsx,y) - 2.*sina.*cosb.*dsyxdspq.*dsyxdsij);
end

end

function dfyxpq = SRVFA_dfyxpq(sina,sinb,cosa,cosb,sigsy,sigsx,y)
tmp = (sigsy./sigsx.*sina.*cosb - sinb.*cosa);
tmp = reshape(tmp.*tmp,[nR,1,1,nTE]);
dfyxpq = SRVFA_dSyxdSpq(sigsy,sigsx,y).*sina.*sinb.*(cosb-cosa)./tmp;
end

function dfyxij = SRVFA_dfyxij(sina,sinb,cosa,cosb,sigsy,sigsx,y)
tmp = (sigsy./sigsx.*sina.*cosb - sinb.*cosa);
tmp = reshape(tmp.*tmp,[nR,1,1,nTE]);
dfyxij = SRVFA_dSyxdSij(sigsy,sigsx,y).*sina.*sinb.*(cosb-cosa)./tmp;
end

function dSyxdSij = SRVFA_dSyxdSij(sigsy,sigsx,y) % make this [nR, 3*nTE, 1,nTE]
% [nR,~,~,nTE] = size((sigsy));
sigsx = reshape(sigsx,[nR,1,nTE]);
sigsy = reshape(sigsy,[nR,1,nTE]);
unusedVal = round(sum(1:3) - y - 1);
% need to copy similar idea but only one derivative. figure out how it
% applies
% nonzero when q==k & (p==x | p==y)
masky = zeros(3,1); masky(y) = 1;
maskx = [1;0;0];
maskxx = maskx;
maskyy = masky;
for count = 1:nTE-1
    maskxx = blkdiag(maskxx,maskx);
    maskyy = blkdiag(maskyy,masky);
end
maskxx = ones(nR,1).*reshape(maskxx,[1,3*nTE,nTE]); % [nR, 3*nTE,nTE]
maskyy = ones(nR,1).*reshape(maskyy,[1,3*nTE,nTE]); % [nR, 3*nTE,nTE]
tmpx = -sigsy./(sigsx.*sigsx).*maskxx;
tmpy = 1./sigsx.*maskyy;
if sum(tmpx.*tmpy,'all') ~=0
    error('Check matrices');
end
dSyxdSij = reshape(tmpx + tmpy,[nR,3*nTE,1,nTE]);

end

function dSyxdSpq = SRVFA_dSyxdSpq(sigsy,sigsx,y)
% [nR,~,~,nTE] = size((sigsy));
tmp = SRVFA_dSyxdSij(sigsy,sigsx,y);
dSyxdSpq = reshape(tmp,[nR,1,3*nTE,nTE]);
end

function dSyxdSpqdSij = SRVFA_dSyxdSpqdSij(sigsy,sigsx,y)
% x always is 1
% [nR,~,~,nTE] = size((sigsy));
sigsx = reshape(sigsx,[nR,1,1,nTE]);
sigsy = reshape(sigsy,[nR,1,1,nTE]);
% if q ~= k, then zero or if (p ~= x and p ~= y) or if q ~=j
% if above and p == y, then i == x (1) or it is zero
% so it is nonzero only in the following cases:
%   q == k == j
% AND
%   for p == y
%     i == x
%   for p == x
%     i == x or i == y
% so make a mask for those cases:
unusedVal = round(sum(1:3) - y - 1);
% maskkqj = ones(3,3);
% maskkqj(:,unusedVal) = 0; maskkqj(unusedVal,:) = 0;
% maskkqj(y,y) = 0; % based on rules above
maskkqj = [1 0 0; 0 0 0; 0 0 0];
masktmp2 = maskkqj;
masktmp = ones(3,1);
for count = 1:(nTE-1)
    maskkqj = blkdiag(maskkqj,masktmp2);
    masktmp = blkdiag(masktmp,ones(3,1));
end
maskkqj1 = ones(nR,1).*reshape(maskkqj.*reshape(masktmp,[3*nTE,1,nTE]),[1,3*nTE,3*nTE,nTE]);
if y == 3
    % do circshift 2 
    maskkqj2 = circshift(maskkqj1,2,2);
    maskkqj3 = circshift(maskkqj1,2,3);
elseif y == 2
    % do circshift 1
    maskkqj2 = circshift(maskkqj1,1,2);
    maskkqj3 = circshift(maskkqj1,1,3);
else 
    error('y must be 2 or 3');
end

tmp1 = 2.*sigsy./(sigsx.*sigsx.*sigsx) .* maskkqj1;
tmp2 = -1./(sigsx.*sigsx);
tmp3 = tmp2 .*maskkqj3;
tmp2 = tmp2.*maskkqj2;

% if sum(tmp1.*tmp2.*tmp3,'all') ~= 0
%     error('Check matrices');
% end
% now I can just do the calculation and apply the mask

dSyxdSpqdSij = tmp1 + tmp2 + tmp3;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dedSijSpq = SRVFA_dedSijSpq(TEv,T2star,T2stard,dT2stard_dSij,dT2star_dSij,sigs1,sigs2,sigs3)
% whole size is [nR, 3*nTE (ij), 3*nTE (pq), nTE (k)]
% nTE = length(TEv);
% nR = length(T2star);
TEv = reshape(TEv,[1,1,1,nTE]);
% T2star is [nR x 1] vector
dT2stard_dSpq = reshape(squeeze(dT2stard_dSij),[nR,1,3*nTE,1]);
dT2star_dSpq = reshape(squeeze(dT2star_dSij),[nR,1,3*nTE,1]);

term1 = (TEv.*TEv).*exp(-TEv.*(1./T2stard - 1./T2star)).*(dT2stard_dSij./(T2stard.*T2stard) - dT2star_dSij./(T2star.*T2star));
term1 = term1.* (1./(T2stard.*T2stard).*dT2stard_dSpq - 1./(T2star.*T2star).*dT2star_dSpq);

dT2stard_dSpqdSij = SRVFA_dT2stard_dSpqdSij(TEv,sigs3);
dT2star_dSpqdSij = SRVFA_dT2star_dSpqdSij(TEv,sigs1,sigs2);
tmp = -2.*dT2stard_dSij./(T2stard.*T2stard.*T2stard).*dT2stard_dSpq + dT2stard_dSpqdSij./(T2stard.*T2stard) + ...
    2.*dT2star_dSij./(T2star.*T2star.*T2star).*dT2star_dSpq - dT2star_dSpqdSij./(T2star.*T2star);

term2 = tmp.*(TEv.*exp(-TEv.*(1./T2stard-1./T2star)));

dedSijSpq = term1+term2;
end

function dT2star_dSij = SRVFA_dT2star_dSij(TEv,sigs1,sigs2)
% [nR,nTE] = size(sigs1);
sigs3use = zeros(size(sigs1)).*NaN;
sigsuse = cat(2,reshape(sigs1,nR,1,nTE),reshape(sigs2,nR,1,nTE),reshape(sigs3use,nR,1,nTE)); % sigs = [nR,3,nTE]
TEv = reshape(TEv,1,1,nTE);
lnsigs = log(sigsuse);

% sumTE = sum(TEv,'all');
top = (nTE.*TEv - sumTE)./sigsuse .*(nTE.*sum(TEv.*TEv,'all') - sumTE.*sumTE);
bot = (nTE.*sum(TEv.*lnsigs,3) - sumTE.*sum(lnsigs,3));
bot = bot.*bot;
dT2star_dSij = .5.*top./bot;
dT2star_dSij(:,3,:) = 0;
dT2star_dSij = reshape(dT2star_dSij,[nR,3*nTE]);
end

function dT2stard_dSij = SRVFA_dT2stard_dSij(TEv,sigs3)
% [nR,nTE] = size(sigs3);
sss1 = zeros(nR,nTE).*NaN;
sss2 = sss1;
sigs = cat(2,reshape(sss1,nR,1,nTE),reshape(sss2,nR,1,nTE),reshape(sigs3,nR,1,nTE)); % sigs = [nR,3,nTE]
TEv = reshape(TEv,1,1,nTE);
lnsigs = log(sigs);

% sumTE = sum(TEv,'all');
top = (nTE.*TEv - sumTE)./sigs .*(nTE.*sum(TEv.*TEv,'all') - (sumTE.*sumTE));
bot = (nTE.*sum(TEv.*lnsigs,3) - sumTE.*sum(lnsigs,3));
bot = bot.*bot;
dT2stard_dSij = .5.*top./bot;
dT2stard_dSij(:,1:2,:) = 0;
dT2stard_dSij = reshape(dT2stard_dSij,[nR,3*nTE]);
end

function dT2stard_dSpqdSij = SRVFA_dT2stard_dSpqdSij(TEv,sigs3)
% only nonzero when i == p
% [nR,~,~,nTE] = size(sigs3);
sss1 = zeros(size(sigs3));
sss2 = sss1;
sigspq = cat(3,reshape(sss1,nR,1,1,nTE),reshape(sss2,nR,1,1,nTE),reshape(sigs3,nR,1,1,nTE)); % sigs = [nR,1,3,nTE]

TEv = reshape(TEv,[1,1,1,nTE]);
TEq = reshape(repelem(TEv(:),3),[1,1,3*nTE,1]); % along the pq dimension
TEj = reshape(TEq,[1,3*nTE,1,1]);
lnsigs = log(sigspq);

% sumTE = sum(TEv,'all');
tmp1 = ((sumTE.*sumTE) - nTE.*sum(TEv.*TEv,'all'))./(nTE.*sum(TEv.*lnsigs,4) - sumTE.*sum(lnsigs,4)).^2; % size = [nR, 1, 3, 1]
tmp1 = repmat(tmp1,[1,1,1,nTE]);
tmp1 = reshape(tmp1,[nR,1,3*nTE,1]); % size = [nR,1,3*nTE,1]

tmp2 = (nTE.*TEq - sumTE)./((reshape(sigspq.*sigspq,[nR,1,3*nTE])));
masktmp = SRVFA_dSpqdSij(nTE);
tmp2 = tmp2.*masktmp;
sigspq2 = reshape(sigspq,[nR,3,nTE]);
Spj = repmat(repelem(permute(sigspq2,[1 3 2]),1,3,1),[1 1 nTE]);
tmp3top = (2./Spj.*SRVFA_dSpjdSij(nTE).*(nTE.*TEj - sumTE)).*(nTE.*TEq - sumTE)./reshape(sigspq,[nR,1,3*nTE]);
tmp3bot = nTE.*repmat(sum(TEv.*lnsigs,4),1,1,nTE,1) - sumTE.*repmat(sum(lnsigs,4),1,1,nTE,1);

dT2stard_dSpqdSij = tmp1.*(tmp2 + tmp3top./tmp3bot);
% dT2stard_dSpqdSij
dT2stard_dSpqdSij(:,:,1:3:end,:) = 0;
dT2stard_dSpqdSij(:,:,2:3:end,:) = 0;

end

function dT2star_dSpqdSij = SRVFA_dT2star_dSpqdSij(TEv,sigs1,sigs2)
% only nonzero when i == p
% [nR,~,~,nTE] = size(sigs1);
sss3 = zeros(size(sigs1));
sigspq = cat(3,reshape(sigs1,nR,1,1,nTE),reshape(sigs2,nR,1,1,nTE),reshape(sss3,nR,1,1,nTE)); % sigs = [nR,1,3,nTE]

TEv = reshape(TEv,[1,1,1,nTE]);
TEq = reshape(repelem(TEv(:),3),[1,1,3*nTE,1]); % along the pq dimension
TEj = reshape(TEq,[1,3*nTE,1,1]);
lnsigs = log(sigspq);

% sumTE = sum(TEv,'all');
tmp1 = ((sumTE.*sumTE) - nTE.*sum(TEv.*TEv,'all'))./(nTE.*sum(TEv.*lnsigs,4) - sumTE.*sum(lnsigs,4)).^2; % size = [nR, 1, 3, 1]
tmp1 = repmat(tmp1,[1,1,1,nTE]);
tmp1 = reshape(tmp1,[nR,1,3*nTE,1]); % size = [nR,1,3*nTE,1]
tmp2 = (nTE.*TEq - sumTE)./((reshape(sigspq.*sigspq,[nR,1,3*nTE]))) .*SRVFA_dSpqdSij(nTE);
sigspq2 = reshape(sigspq,[nR,3,nTE]);
Spj = repmat(repelem(permute(sigspq2,[1 3 2]),1,3,1),[1 1 nTE]);
tmp3top = (2./Spj.*SRVFA_dSpjdSij(nTE).*(nTE.*TEj - sumTE)).*(nTE.*TEq - sumTE)./reshape(sigspq,[nR,1,3*nTE]);
tmp3bot = nTE.*repmat(sum(TEv.*lnsigs,4),1,1,nTE,1) - sumTE.*repmat(sum(lnsigs,4),1,1,nTE,1);

dT2star_dSpqdSij = 0.5.*tmp1.*(tmp2 + tmp3top./tmp3bot);
dT2star_dSpqdSij(:,:,3:3:end,:) = 0;
end

function dSpqdSij = SRVFA_dSpqdSij(nTE)
% basically, if the signals are equal, then the answer is 1, if not, then
% is 0. This is just a mask basically
TE3 = 3*nTE;
dSpqdSij = reshape(eye(TE3),[1,TE3,TE3]);
end

function dSpjdSij = SRVFA_dSpjdSij(nTE)
% baslically this is another mask, but a little different
% when i == p, it is 1, else it is 0. so it is uniform across q dimension
TE3 = 3*nTE;
dSpjdSij = reshape(repmat(eye(3),nTE,nTE),[1,TE3,TE3]);
end
end