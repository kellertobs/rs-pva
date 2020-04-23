% varimax:  rs-pva subroutine to perform varimax rotation on k endmember
%           compositions and corresponding mixing proportions.
%
%   [Appvm] = varimax(App)
% 
%   Performs varimax rotation: reverses axes 180 degrees to maximize positive 
%   mixing proportions. The procedure follows the algorithm as spelled out in
%   Harman (1960) in Chapter 14, section 4. The code is modified from varimax.m 
%   included in MATLAB Appendix of Reyment and Joreskog (1993). The notation 
%   follows Harman (1960).
%
%   App    : input mixing proportions in transformed factor space reduced to k endmembers
%   Appvm  : output varimax rotated mixing proportions, maximised positive entries
%            and row sums normalised to unity.
%
% created  : 2020-03-21  Tobias Keller, University of Glasgow
% based on : 1994-09-12  Glenn Johnson, University of Utah
% license  : GNU General Public License v3.0


function [Appvm] = varimax(App)

% get matrix size and scale
[m,n] = size(App);
hj    = sqrt(diag(App*App'));

% compute variances of squared mixing proportions
Ah = App./hj;
V0 = m*sum(sum(Ah.^4))-sum(sum(Ah.^2).^2);

res = 1;
tol = 1e-4;
it  = 1;
while res >= tol % iterate until residual below tolerance
    for i=1:n-1  % cycles through 2 factors at a time.
        jl=i+1;
        for j=jl:n
            xj  = App(:,i)./hj;   % notation here closely
            yj  = App(:,j)./hj;   % follows harman
            
            uj  = xj.*xj-yj.*yj;
            vj  = 2*xj.*yj;
            
            a   = sum(uj);
            b   = sum(vj);
            c   = uj'*uj-vj'*vj;
            d   = 2*uj'*vj;
            
            num = d-2*a*b/m;
            den = c-(a^2-b^2)/m;
            
            phi   = atan2(num,den)/4;
            
            if abs(phi) > tol
                Xj  = cos(phi)*xj+sin(phi)*yj;
                Yj  =-sin(phi)*xj+cos(phi)*yj;
                
                App(:,i) = Xj.*hj;
                App(:,j) = Yj.*hj;
            end
        end
    end
    
    Appvm = App;
    Ah    = App./hj;
    V     = m*sum(sum(Ah.^4))-sum(sum(Ah.^2).^2);
    res   = (V-V0)./V0;
    V0    = V;
    it    = it+1;
    if it > 20; error('Varimax rotation did not converge after 20 iterations!'); end
end

% End varimax

% Maximize positive mixing proportions
for j=1:n
    if abs(min(Appvm(:,j))) > abs(max(Appvm(:,j)))
        Appvm(:,j) = (-1)*Appvm(:,j);
    end
end

% Normalize varimax mixing proportions to unity row sum
Appvm  = Appvm./sum(Appvm,2);
