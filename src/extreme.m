% extreme:  rs-pva subroutine to find k mutually extreme samples in dataset
%   
%   [App0,Fpp0] = extreme(Appvm,Fppvm,Xns,DGN)
%
%   Function determines the k mutually extreme samples in dataset in row-sum 
%   normalised data Xns. Input Appvm and Fppvm are arrays with scaled and 
%   varimax rotated mixing proportions and principal components for a mixing
%   model with k endmembers. The routine passes back scaled initial guesses 
%   for mixing proportions App0 and k endmember compositions Fpp0 as input
%   for the rs-pva optimisation routine. Notation follows Miesch (1976).
%
%   Appvm  : input varimax rotated mixing proportions in transformed factor space
%   Fppvm  : input varimax rotated endmember compositions in transformed factor space
%   Xns    : row-sum normalised data in measurment space
%   DGN    : input rs-pva diagnostics structure for data format metrics
%   App0   : input starting guess mixing proportions in transformed factor space
%   Fpp0   : input starting guess endmember compositions in transformed factor space

% created  : 2020-03-21  Tobias Keller, University of Glasgow
% based on : 1994-09-12  Glenn Johnson, University of Utah
% license  : GNU General Public License v3.0


function    [App0,Fpp0] = extreme(Appvm,Fppvm,Xns,DGN)

k    = DGN.k;
rows = zeros(k,1);

% Get k unique samples with extreme mixing proportions from A''
[~,Irs] = sort(abs(Appvm),1,'descend','ComparisonMethod','abs');
for i = 1:k
    j = 1;
    if any(Irs(j,i)==Irs(j,(1:k)<i))
        rows(i) = Irs(2,i);
    else
        rows(i) = Irs(1,i);
    end
end

% Calculate oblique transform matrix O0 (rows of A'' for k unique max. mixing proportions)
O0 = Appvm(rows,:);

% Calculate initial oblique endmembers and mixing proportions
Fpp0 = O0*Fppvm;
App0 = Appvm/O0;

% Save initial oblique transform matrix in case EXTREME algorithm does not converge
O0i    = O0;
rowsin = rows;

% Begin EXTREME iterations
maxit   = 20;
rowsold = 0.*rows;
iter    = 1;
while any(rows ~= rowsold)
    
    rowsold = rows;

    % scale matrices to measurement units
    [A0,~]  = scaleup(Fpp0,App0,Xns);
    
    % Get k unique samples with extreme mixing proportions from A0
    [~,Irs] = sort(A0,1,'descend');
    for i = 1:k
        j = 1;
        if any(Irs(j,i)==Irs(j,(1:k)<i))
            rows(i) = Irs(2,i);
        else
            rows(i) = Irs(1,i);
        end
    end
    
    % Extract updated oblique transform matrix
    O0 = Appvm(rows,:);
    
    % calculate updated endmembers and proportions
    Fpp0 = O0*Fppvm;
    App0 = Appvm/O0;
    
    % Report convergance
    if all(rowsold == rows)
        disp(['  - Extreme samples determination converged after ',int2str(iter),' iterations:'])
        disp(['    Initial polytope vertices will use samples: ',num2str(rows.')]);
        disp(' ');
    end
    
    % Stop iterating and use initial extreme mix. proportions when maxits reached
    if iter == maxit
        disp(['  ! Extreme samples determination did not converge after ',int2str(iter),' iterations:'])
        disp(['    Initial polytope vertices will use samples: ',num2str(rowsin.')]);
        disp(' ');
        O0      = O0i;
        Fpp0    = O0*Fppvm;
        App0    = Appvm/O0;
        break;
    end
    
    iter = iter+1;
end

end  % end of function
