
function [Abf,Fbf,Xbf,DGN] = solve(A0,F0,X0,DGN,VNAME)

% initialise variables and parameters
Xbf     = X0;
Xf      = X0;
Af      = A0;
Ff      = F0;
Ir      = DGN.Ir;
Ii      = DGN.Ii;
rm      = DGN.rm;
As      = sort(Af);
misfit  = 1;
bestfit = 1;
rmvcrt  = 0;
Abf     = Af;
Fbf     = Ff;
it      = 1;
rng(15);

% adjust convergence tolerances and step parameters
dft = [1e-3, 0.05, 1000, 0.5, 0.25];
par = input(['->  Adjust convergence parameters as list: [msftol,negtol,rmvtol,alpha,beta] \n' ...
             '    msftol: misfit tolerance for convergence          (dft = ',num2str(dft(1)),') \n' ...
             '    negtol: tolerance for negative mixing proportions (dft = ',num2str(dft(2)),') \n' ...
             '    rmvcnt: unsuccesful iterations to outlier removal (dft = ',num2str(dft(3)),') \n' ...
             '    alpha : step size for directed iterative update   (dft = ',num2str(dft(4)),') \n' ...
             '    beta  : step size for random iterative update     (dft = ',num2str(dft(5)),') \n']);
if isempty(par); par = dft; end

msftol  = par(1);
negtol  = par(2);
rmvcnt  = par(3);
alpha   = par(4);
beta    = par(5);

% set preference for plotting progress
dft = 1;
plt = input('->  Do you wish to plot progress of optimisation (yes = 1, dft) or not (no = 0)? \n');
if isempty(plt); plt = dft; end

% iterative optimisation algorithm, minimise misfit function below tolerance
while misfit > msftol
    
    % update mixing proportions (no smaller than tolerance, unity sum)
    D  =   - alpha .* min(0,As(1,:)+negtol);  % directed increment towards less negative mix prop
    D  = D + beta  .* (rand(1,DGN.k)-0.5).*min(1,bestfit);  % random increment
    Af(Ii,:) = (Abf(Ii,:) + D)./(1+sum(D));  % update mixing proportions
    As = sort(Af);  % get sorted mix prop
    
    % update EM compositions (no negatives, 100% sum)
    Ff = (Af(Ii,:).'*Af(Ii,:))\(Af(Ii,:).'*X0(Ii,:));  % Lsq-solve for EM comp
    Ff = max(0,Ff);  % zero out negative EM comp
    Ff = Ff./sum(Ff,2).*100;  % normalise EM comp to 100% sum
    
    % update data fit
    Xf(Ii,:) = Af(Ii,:)*Ff;  
    
    % calculate misfit
    msft_dtft = norm((Xf-X0)./mean(X0),2);  % calculate data misfit 
    msft_ngpr = norm(min(0,Af(Ii,:)+negtol),2);  % calculate neg. mix. prop. misfit
    misfit    = msft_dtft + msft_ngpr;  % get combined misfit
    
    % update best fit and report progress
    if misfit < bestfit
        
        if plt % visualise best fit k-EM model and fitted data
            DGN.Ii = Ii; DGN.Ir = Ir; DGN.rm = rm;
            visualise({X0,Xf,Fbf},{'data','fitted data','fitted EM'},['it ',num2str(it),';  misfit = ',num2str(bestfit,4),';  removed ',num2str(length(Ir)),' outliers;'],DGN,VNAME)
        end
        
        % report new best fit and diagnostics
        if                  it <  10
            disp(['    ---  it    ',int2str(it),';  misfit ',num2str(bestfit,'%1.3e'),';  removed ',int2str(rm),' outliers;'])
        elseif it >= 10  && it < 100
            disp(['    ---  it   ',int2str(it),';  misfit ',num2str(bestfit,'%1.3e'),';   removed ',int2str(rm),' outliers;'])
        elseif it >= 100 && it < 1000
            disp(['    ---  it  ',int2str(it),';  misfit ',num2str(bestfit,'%1.3e'),';  removed ',int2str(rm),' outliers;'])
        else
            disp(['    ---  it ',int2str(it),';  misfit ',num2str(bestfit,'%1.3e'),';  removed ',int2str(rm),' outliers;'])
        end
        
        % update best fit solution
        bestfit = misfit;
        Fbf     = Ff;
        Abf     = Af;
        Xbf     = Xf;
        
        % reset outlier removal criterion
        rmvcrt = 0;
    else
        % increment outlier removal criterion
        rmvcrt = rmvcrt+1;
    end
    
    % remove data outlier if criterion reached
    if rmvcrt > rmvcnt
        % determine outlier                         % add tiny random number to avoid zero criterion
        minA      = min(min(0,Abf+negtol),[],2)     + rand(length(Abf),1).*1e-16;  % row-wise lowest neg. mix. prop.
        maxR      = max(abs(Xbf-X0)./mean(X0),[],2) + rand(length(Xbf),1).*1e-16;  % row-wise highest data misfit
        ir        = find(minA==min(minA)).*(abs(min(minA))>=max(maxR)) ... % select lowest negative mix. prop.
                  + find(maxR==max(maxR)).*(abs(min(minA))< max(maxR));    % or highest data misfit
        
        % remove outlier and update removal diagnostics
        Ir         = [Ir;ir];
        Ii(Ii==ir) = [];
        rm         = length(Ir);
        Af(Ir,:)   = 0;
        Abf(Ir,:)  = 0;
        Xf(Ir,:)   = X0(Ir,:);
        Xbf(Ir,:)  = X0(Ir,:);
        As         = sort(Af);  % get sorted mix prop

        % stop optimisation if >10% data identified as outliers
        if rm > DGN.m/10; error('Optimisation failed: too many samples identified as outliers.'); end
        
        % reset outlier removal criterion
        rmvcrt    = 0;
    end
    
    it = it+1; % increment iteration count
end

% bundle diagnostics for output
DGN.it = it;
DGN.bf = bestfit;
DGN.rm = rm;
DGN.Ir = Ir;
DGN.Ii = Ii;

end