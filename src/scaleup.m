% scaleup:  rs-pva subroutine to scale endmember compositions and mixing 
%           proportions from transformed factor to measurement space.
%
%   [A,F] = scaleup(Fpp,App,Xns)
% 
%   Fpp    : input endmember compositions in transformed factor space
%   App    : input mixing proportions in transformed factor space
%   Xns    : input row-sum normalised data in measurment space
%   F      : output endmember compositions scaled up to measurement units (row sum = 100%)
%   A      : output mixing proportions scaled up to measurement units (row sum = 1)
%
% created  : 2020-03-21  Tobias Keller, University of Glasgow
% based on : 1994-09-12  Glenn Johnson, University of Utah
% license  : GNU General Public License v3.0


function    [A,F] = scaleup(Fpp,App,Xns)

% calculate scale factors
SUMxmin = sum(min(Xns)); % numerator of Miesch equation 9
SUMfran = sum(Fpp.*(max(Xns)-min(Xns)),2); % denominator of Miesch equation 9
skpc    = (100-SUMxmin)./SUMfran; % Miesch equation (9)

% calculate intermediate scaled oblique endmembers [F']
Fp = Fpp.*skpc;

% calculate unscaled oblique endmembers in measurement units [F]
F  = Fp.*(max(Xns)-min(Xns)) + min(Xns);

% calculate intermediate scaled mixing proportions [A']
Ap = App./skpc.';

% calculate unscaled mixing proportions in measurement units [A]
A  = Ap./sum(Ap,2);

