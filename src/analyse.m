% analyse:  rs-pva subroutine to analyse dataset to determine optimal number 
%           of endmembers for mixing model
%
%   [Xns,Xnr,Xsq,F,A,Fpp,App,DGN] = analyse(X)
%
%   X      : input raw data array in measurement units
%   Xns    : output data array transformed to 100% row sum
%   Xns    : output data array transformed to unit range
%   Xns    : output data array transformed by square-root of sum of squares
%   F      : output endmember compositions scaled up to measurement units (row sum = 100%)
%   A      : output mixing proportions scaled up to measurement units (row sum = 1)
%   Fpp    : output endmember compositions in transformed factor space
%   App    : output mixing proportions in transformed factor space
%   DGN    : input/output rs-pva diagnostics structure for various data format and statistical metrics
%
% created  : 2020-03-21  Tobias Keller, University of Glasgow
% based on : 1994-09-12  Glenn Johnson, University of Utah
% license  : GNU General Public License v3.0


function [Xns,Xnr,Xsq,F,A,Fpp,App,DGN] = analyse(X,DGN)

% Try alternative ways of calculating CDtab to get it to match
% SAWVEC CD table.

% get m: # samples, n: # variables
m = DGN.m;
n = DGN.n;

% Calculate descriptive statistics of untransformed data
untrstat = ones(n,5);
untrstat(:,1) = 1:n;
untrstat(:,2) = mean(X);
untrstat(:,3) = std (X);
untrstat(:,4) = min (X);
untrstat(:,5) = max (X);

% Transform data to constant row sum
Xns = X./sum(X,2).*100;

% Calculate descriptive statistics of constant row sum transformed data
nsumstat = ones(n,5);
nsumstat(:,1) = 1:n;
nsumstat(:,2) = mean(Xns);
nsumstat(:,3) = std (Xns);
nsumstat(:,4) = min (Xns);
nsumstat(:,5) = max (Xns);

% Transform data to normalised range [0,1]
Xnr = (Xns-min(Xns))./(max(Xns)-min(Xns));

% Calculate descriptive statistics of normalised range transformed data
nrngstat = ones(n,5);
nrngstat(:,1) = 1:n;
nrngstat(:,2) = mean(Xnr);
nrngstat(:,3) = std (Xnr);
nrngstat(:,4) = min (Xnr);
nrngstat(:,5) = max (Xnr);

% transform data to row-wise normalised by square root of sum of squares
Xsq = Xnr./(sum(Xnr.^2,2).^0.5);

% perform singular value decomposition and calculate variance diagnostics: Q = u*s*v'
[u,s,v] = svd(Xsq,0);
sv      = ones(n,4);
sv(:,1) = 1:n;
sv(:,2) = diag(s).^2;
sv(:,3) = diag(s).^2./sum(s(:).^2).*100;
sv(:,4) = cumsum(sv(:,3));

% CALCULATE MIESCH SCALED PRINCIPAL COMPONENTS [F'', A'']

% Put matrices in Miesch (1976) terminology: Q = App*Fpp
Fpp = v.';  % scaled principal endmember components
App = u*s;  % scaled principal mixing proportions

% Scale back up to measurement units
[A,F] = scaleup(Fpp,App,Xns);

% % Calculate scale factor
% SUMxmin = sum(min(Xns)); % numerator of Miesch equation 9
% SUMfran = sum(Fpp.*(max(Xns)-min(Xns)),2); % denominator of Miesch equation 9
% skpc    = (100-SUMxmin)./SUMfran; % Miesch equation (9)
% 
% % Calculate [F'] in Miesch equation 12
% Fp = Fpp.*skpc;
% 
% % Calculate [F] in Miesch equation 16
% F = Fp.*(max(Xns)-min(Xns)) + min(Xns);
% 
% % Calculate Miesch [A']
% Ap = App./skpc.';
% 
% % Calculate Miesch [A]
% A  = Ap./sum(Ap,2);

% Calculate Goodness of Fit Statistics: Coefficients of Determination
CDtab = ones(n,n-1);
for z = 1:n-1
    kk   = z+1;
    Appr = App(:,1:z+1); % Reduced dimensional App
    Fppr = Fpp(1:z+1,:); % Reduced dimensional Fpp
    Xsqr = Appr*Fppr;
    
    disp(['    Calculating Varimax Matrices for ',int2str(z+1),' Components']);
    Appvm = varimax(Appr); % get varimax rotated, scaled mixing proportions
    
    Fppvm   = (Appvm'*Appvm)\(Appvm'*Xsqr);  % get scaled endmember compositions
    [Appr,Fppr] = scaleup(Fppvm,Appvm,Xns);      % scale back endmembers and proportions to measurement units
    Xnsf    = Appr*Fppr;                         % calculate row-sum normalised fitted data
    
    resid = Xns-Xnsf;  % get residual between model fit and data
    for i = 1:n
        S = std(Xns(:,i))^2;
        D = std(resid(:,i))^2;
        CDtab(i,z) = 1-D/S;  % calculate misfit
    end
end
disp(' ');

% Calculate Communalities
COMMtab = zeros(m,n-1);
for j = 1:n-1
    Appr = App(:,1:j+1);
    COMMtab(:,j) = diag(Appr*Appr');
end

% Assemble diagnostics structure
DGN.stat.untr = untrstat;
DGN.stat.nsum = nsumstat;
DGN.stat.nrng = nrngstat;
DGN.SVD       = sv;
DGN.CDvar     = CDtab;
DGN.CMsmp     = COMMtab;
DGN.m         = m;
DGN.n         = n;
