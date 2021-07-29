function [refAng,dynAng] = findBestSRVFAAngle(TR,T1,T2star,T1dovT1,zeta,...
    TE,aveDims,doT2StarCorrection,doT2StarAll,numBaselineAvg,numDynamicAvg)
% Description: findBestSRVFAAngle - finds the "optimal" flip angle combination to
%   minimize the relative standard deviation of T1 maps produced by the
%   single reference variable flip angle method.
%
% Inputs:
%     TR (length nTR vector): list of repetition times
%     T1 (length nT1 vector): list of T1,baseline values
%     T2star (length nT2s vector): list of T2star values
%     T1dovT1 (length nT1d vector): list of ratios of dynamic T1 / baseline
%         T1 to include in the minimization. For example, with a starting T1 or
%         200, a T1dovT1 ratio of 1:0.2:2 would include dynamic T1s of [200,
%         240, 280, 320, 360, 400]
%     zeta (length nZeta vector): list of sensitivity ratios (zeta)
%     TE (length nTE vector): list of echo times in same units as T2star.
%     aveDims (string vector): dims over which to average the input parameter
%         sets. Valid values are as follows. Double quotes must be used for
%         this to be interpreted properly:
%         ["TR","T1","T1dovT1","TE","T2star","zeta"]
%     doT2StarCorrection (0 or 1): 
%         0 = do SR-VFA method
%         1 = do SR-VFA-T2* method
%     doT2StarAll (0 or 1):
%         0 = do T2star cor. without echo combination
%         1 = do T2star cor. with inverse variance weighted echo combination
%     numBaselineAvg (integer): number of baseline scan averages
%     numDynamicAvg (integer): number of dynamic scan averages 
%
% Outputs:
%     refAng (nTR x nT1 x nT2s x nT1d x nZeta x nTE matrix): reference flip
%         angle for minimization of T1d variance
%     dynAng (nTR x nT1 x nT2s x nT1d x nZeta x nTE matrix): dynamic flip
%         angle for minimization of T1d variance
%
%     Note: The length of each averaged dimension will be 1
%
% Author: 
%     Michael Malmberg
%     University of Utah
%     2021
%
% For an original description of the SR-VFA method, see 
% 
% Svedin, B. T., Payne, A., & Parker, D. L. (2019). 
% Simultaneous proton resonance frequency shift thermometry and T1 
% measurements using a single reference variable flip angle T1 method. 
% Magnetic Resonance in Medicine, 81(5), 3138�3152. 
% https://doi.org/10.1002/mrm.27643
%
% For the T2* correction and a detailed analysis on the variance of the
% above method (from whence the equations used herein come), see
%
% Malmberg, M. Odeen, H., & Parker, D. L. (2021).
% Effects of T2* on accuracy and precision of single-reference variable
% flip angle T1-mapping.
%

%% Do error checking
if ~exist('numBaselineAvg','var') || isempty(numBaselineAvg)
    numBaselineAvg = 1;
end
if ~exist('numDynamicAvg','var') || isempty(numDynamicAvg)
    numDynamicAvg = 1;
end
baselineAvg = numBaselineAvg;
dynamicAvg = numDynamicAvg;
try 
    if mod(baselineAvg,1) ~=0 || mod(dynamicAvg,1) ~=0 || baselineAvg < 1 || dynamicAvg < 1
        error('baseline and dynamic averages must be positive integers');
    end
catch
    error('baseline and dynamic averages must be positive integers');
end

if ~exist('doT2StarCorrection','var') || isempty(doT2StarCorrection)
    doT2StarCorrection = false;
end
if ~exist('doT2StarAll','var') || isempty(doT2StarCorrection) || doT2StarCorrection == false
    doT2StarAll = false;
end
if ~exist('aveDims','var') 
    aveDims = [];
end
valids = ["TR","T1","T1dovT1","TE","T2star","zeta"];
aveDims = string(aveDims);
for ii = 1:length(aveDims)
    if ~isempty(aveDims) && (isnumeric(aveDims(ii)) || ~ismember(aveDims(ii),valids))
        disp(valids);
        error('aveDims can only include the above listed items');
    end
end

%% Initialize Variables and functions
M0 = 1;
sigma = 1e-2; % arbitrary noise value is chosen, since it (most of the time) doesn't affect the calculation's result.
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','FunctionTolerance',1e-7,'Display','off');
stdT1Find = @(x,xdata) (SRVFA_SimStd(M0,xdata{1},xdata{2},xdata{3},x(1),x(2),xdata{4},xdata{5},xdata{6},sigma,aveDims,doT2StarCorrection,doT2StarAll,baselineAvg,dynamicAvg));
stdT1Find2 = @(x,xdata) (SRVFA_SimStd(M0,xdata{1},xdata{2},xdata{3},x(1),x(2),xdata{4},xdata{5},xdata{6},sigma*100,aveDims,doT2StarCorrection,doT2StarAll,baselineAvg,dynamicAvg));
stdT1Find3 = @(x,xdata) (SRVFA_SimStd(M0,xdata{1},xdata{2},xdata{3},x(1),x(2),xdata{4},xdata{5},xdata{6},sigma*10,aveDims,doT2StarCorrection,doT2StarAll,baselineAvg,dynamicAvg));
stdT1Find4 = @(x,xdata) (SRVFA_SimStd(M0,xdata{1},xdata{2},xdata{3},x(1),x(2),xdata{4},xdata{5},xdata{6},sigma*1000,aveDims,doT2StarCorrection,doT2StarAll,baselineAvg,dynamicAvg));
% xdata is as follows:
%   TR,T1,T1d,TE,T2s,T2sd
% x is as follows:
%   FA1, FA2
nTR = length(TR);
nT1 = length(T1);
nT2s = length(T2star);
nT1d = length(T1dovT1);
nZeta = length(zeta);
nTE = length(TE);

ydata = 0;

if doT2StarCorrection && doT2StarAll && ~ismember("TE",aveDims)
    aveDims = [aveDims, "TE"];
end
for ii = 1:length(aveDims)
    switch aveDims(ii)
        case "TR"
            nTR = 1;
            xdata{1} = TR;
        case "T1"
            nT1 = 1;
            xdata{2} = T1;
        case "T2star"
            nT2s = 1;
            xdata{3} = T2star;
        case "T1dovT1"
            nT1d = 1;
            xdata{4} = T1dovT1;
        case "zeta"
            nZeta = 1;
            xdata{5} = zeta;
        case "TE"
            nTE = 1;
            xdata{6} = TE;
    end
end

if doT2StarCorrection
    xdata{6} = TE;
    if ~ismember('TE',aveDims) && ~doT2StarAll
        ydata = zeros(1,1,1,1,1,1,1,1,length(TE));
    end
end

%% Loop through minimization on all parameter sets
for ii = 1:nTR
    if ~ismember("TR",aveDims)
        xdata{1} = TR(ii);
    end
    for jj = 1:nT1
        if ~ismember("T1",aveDims)
            xdata{2} = T1(jj);
        end
        x0(1) = 0.5.*acosd(exp(-TR(ii)/xdata{2}(end)));
        x0(2) = 2*x0(1)+20;
        for kk = 1:nT2s
            if ~ismember("T2star",aveDims)
                xdata{3} = T2star(kk);
            end
            for ll = 1:nT1d
                if ~ismember("T1dovT1",aveDims)
                    xdata{4} = T1dovT1(ll);
                end
                for mm = 1:nZeta
                    if ~ismember("zeta",aveDims)
                        xdata{5} = zeta(mm);
                    end
                    for nn = 1:nTE
                        if ~ismember("TE",aveDims) && ~(doT2StarCorrection)
                            xdata{6} = TE(nn);
                        end
                        X = lsqcurvefit(stdT1Find,x0,xdata,ydata,[],[],options);
                        % check and make sure that the dynamic flip angle is larger than the reference flip angle
                        if abs(X(1))>abs(X(2)) 
                            X = lsqcurvefit(stdT1Find2,x0,xdata,ydata,[],[],options);
                        end
                        if abs(X(1))>abs(X(2))
                            X = lsqcurvefit(stdT1Find3,x0,xdata,ydata,[],[],options);
                        end
                        if abs(X(1))>abs(X(2))
                            X = lsqcurvefit(stdT1Find4,x0,xdata,ydata,[],[],options);
                        end
                        refAng(ii,jj,kk,ll,mm,nn) = X(1);
                        dynAng(ii,jj,kk,ll,mm,nn) = X(2);
                    end
                end
            end
        end
    end
end

refAng = abs(refAng); dynAng = abs(dynAng);
