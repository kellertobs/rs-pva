function [Appvm,Anormvm] = varimax(A0)

% function [Fvm,Avm]=kv_vmax1(A0)
%
% created : Glenn Johnson 1994-09-12
% modified: Tobias Keller 2020-03-21
%
% Varimax Subroutine
% Performs varimax rotation: reverses axes 180 degrees to maximize positive 
% mixing proportions. The procedure follows the algorithm as spelled out in
% Harman (1960) in Chapter 14, section 4. The code is modified from varimax.m 
% included in MATLAB Appendix of Reyment and Joreskog (1993). The notation 
% follows Harman (1960).

b     = A0;
[m,n] = size(A0);
hjsq  = diag(A0*A0');   % communalities
hj    = sqrt(hjsq);

% compute variances of loadings^2
bh = A0./hj;
V0 = m*sum(sum(bh.^4))-sum(sum(bh.^2).^2);

res = 1;
tol = 1e-4;
it  = 1;
while res >= tol % iterate until residual below tolerance
    for i=1:n-1  % cycles through 2 factors at a time.
        jl=i+1;
        for j=jl:n
            xj  = A0(:,i)./hj;   % notation here closely
            yj  = A0(:,j)./hj;   % follows harman
            
            uj  = xj.*xj-yj.*yj;
            vj  = 2*xj.*yj;
            
            A   = sum(uj);
            B   = sum(vj);
            C   = uj'*uj-vj'*vj;
            D   = 2*uj'*vj;
            
            num = D-2*A*B/m;
            den = C-(A^2-B^2)/m;
            
            phi   = atan2(num,den)/4;
            
            if abs(phi) > tol
                Xj  = cos(phi)*xj+sin(phi)*yj;
                Yj  =-sin(phi)*xj+cos(phi)*yj;
                
                bj1    = Xj.*hj;
                bj2    = Yj.*hj;
                b(:,i) = bj1;
                b(:,j) = bj2;
                
                A0(:,i) = b(:,i);
                A0(:,j) = b(:,j);
            end
        end
    end
    
    Appvm = b;
    bh  = A0./hj;
    V   = m*sum(sum(bh.^4))-sum(sum(bh.^2).^2);
    res = (V-V0)./V0;
    V0  = V;
    it  = it+1;
    if it > 20; error('Varimax rotation did not converge after 20 iterations!'); end
end

% End varimax

% Maximize positive mixing proportions
for j=1:n
    if abs(min(Appvm(:,j))) > abs(max(Appvm(:,j)))
        Appvm(:,j) = (-1)*Appvm(:,j);
    end
end

% Normalized varimax mixing proportions
Anormvm  = Appvm./sum(Appvm,2);
