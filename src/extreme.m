% Script to find k mutually extreme samples in dataset.
% Notation follows Miesch (1976)
%
% created :            Glenn Johnson
% modified: 2020-04-21 Tobias Keller

rowsin = zeros(k,1);
[Y,I]  = sort(abs(Appvm));
kones  = ones(k,1);

% Get samples with max mixing proportions in A'' for each principal component without duplication
for i=1:k
    mm   = m;
    samp = I(mm,i);
    chk  = samp*kones;
    if all(chk-rowsin)
        rowsin(i) = samp;
        chk(i)    = samp+1;
    else
        while ~all(chk-rowsin)
            samp = I(mm-1,i);
            chk  = samp*kones;
            if all(chk-rowsin)
                rowsin(i) = samp;
                chk(i)    = samp+1;
            else
                mm = mm-1;
            end
        end
    end
end

% Calculate oblique transform matrix O0 (rows of A'' for three unique max. mixing proportions)
O0 = Appvm(rowsin,:);

% Calculate initial oblique endmembers and mixing proportions
Fpp0 = O0*Fppvm;
App0 = Appvm/O0;

% Save initial oblique transform matrix in case XRAWC does not converge
O0i     = O0;
rowsini = rowsin;

% Begin XRAWC iterations
niter   = 30;
rowsold = 0.*rowsin;
iter    = 1;
while rowsin ~= rowsold
    
    % Scale Matrices to measurement units
    [A0,~]  = scaleup(Fpp0,App0,Xns);
    rowsold = rowsin;
    
    % Sort mixing proportions
    [Y,I]  = sort(A0);
    rowsin = zeros(k,1);
    kones = ones(k,1);
    
    % Get max loadings samples without duplication
    for i=1:k
        mm   = m;
        samp = I(mm,i);
        chk  = samp*kones;
        if all(chk-rowsin)
            rowsin(i) = samp;
            chk(i)    = samp+1;
        else
            while ~all(chk-rowsin)
                samp = I(mm-1,i);
                chk  = samp*kones;
                if all(chk-rowsin)
                    rowsin(i) = samp;
                    chk(i)    = samp+1;
                else
                    mm = mm-1;
                end
            end
        end
    end
    
    % Extract new oblique transform matrix
    O0 = Appvm(rowsin,:);
    
    % Stop iterating and use initial max proportions when niter reached
    if iter == niter
        disp(['  ! Extreme samples determination did not converge after ',int2str(iter),' iterations:'])
        disp(['    Initial polytope vertices will use maximum varimax proportions: ',num2str(rowsini.')]);
        disp(' ');
        rowsin  = rowsini;
        rowsold = rowsin;
        O0      = O0i;
    else
        disp(['  - Extreme samples determination converged after ',int2str(iter),' iterations:'])
        disp(['    Initial polytope vertices will use samples: ',num2str(rowsin.')]);
        disp(' ');
    end
    
    % calculate updated endmembers and proportions
    Fpp0 = O0*Fppvm;
    App0 = Appvm/O0;
    
    iter = iter+1;
end

[A0,F0] = scaleup(Fpp0,App0,Xns);
