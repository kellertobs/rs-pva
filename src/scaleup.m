function    [A,F] = scaleup(Fpp,App,Xns)
% function [A1,F1]=kvb_scbk(m,n,kk,SUMxmin,ksumstat,Fpp1,App1)
% MATLAB function  kvb_scbk.m
% INPUT MATRICES
%    m = number of samples
%    n = number of variables
%    kk=number of eigenvectors/factors retained
%    SUMxmin  = SUM of the n variable's minimum values in constant row sum matrix
%    ksumstat = constant row sum descriptive statistics matrix
%               which has max and min for each variable.
%    Fpp1 = Scores in transformed factor space
%    App1 = Loadings in transformed factor space
% OUTPUT MATRICES
%    A1  = Scaled Loadings Matrix (EM Mixing Proportions) (Row Sum = 1)
%    F1  = Scaled Scores Matrix (EM Compositions) Rows Sum to 100%

% CALCULATE SCALE FACTORS
SUMxmin = sum(min(Xns)); % numerator of Miesch equation 9
SUMfran = sum(Fpp.*(max(Xns)-min(Xns)),2); % denominator of Miesch equation 9
skpc    = (100-SUMxmin)./SUMfran; % Miesch equation (9)

% Calculate oblique [F']
Fp = Fpp.*skpc;

%Calculate oblique [F]
F = Fp.*(max(Xns)-min(Xns)) + min(Xns);

% Calculate Miesch [A']
Ap = App./skpc.';

% Calculate Miesch [A]
A  = Ap./sum(Ap,2);

