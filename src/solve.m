
function [Abf,Fbf,Xbf,DGN] = solve(A0,F0,X0,DGN,VNAME)

% initialise variables and parameters
[k,n]   = size(F0);
[m,~]   = size(A0);
A0      = A0(DGN.Ii,:);
Xr      = X0(DGN.Ir,:);
X0      = X0(DGN.Ii,:);
Xbf     = X0;
Af      = A0;
Ff      = F0;
Ir      = DGN.Ir;
Ii      = DGN.Ii;
rm      = length(Ir);
As      = sort(Af);
misfit  = 1;
bestfit = 1;
rmvcrt  = 0;
Abf     = Af;
Fbf     = Ff;
it      = 1;

% adjust convergence tolerances and step parameters
dft = [1e-6, 0.05, 1000, 0.5, 0.1];
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
    D  = D + beta  .* (rand(1,k)-0.5).*min(1,bestfit);  % random increment
    Af = (Abf + D)./(1+sum(D));  % update mixing proportions
    As = sort(Af);  % get sorted mix prop
    
    % update EM compositions (no negatives, 100% sum)
    Ff = (Af.'*Af)\(Af.'*X0);  % Lsq-solve for EM comp
    Ff = max(0,Ff);  % zero out negative EM comp
    Ff = Ff./sum(Ff,2).*100;  % normalise EM comp to 100% sum
    
    % update data fit
    Xf = Af*Ff;  
    
    % calculate misfit
    msft_dtft = norm((Xf-X0)./mean(X0),2); % calculate data misfit 
    msft_ngpr = norm(min(0,Af+negtol),2);  % calculate neg. mix. prop. misfit
    misfit    = msft_dtft + msft_ngpr;     % get combined misfit
        
    % update best fit and report progress
    if misfit < bestfit
        
        if plt % visualise best fit k-EM model and fitted data
            FS = {'FontSize',14}; MS = {'MarkerSize',8};
            ii = ceil(sqrt(n-1)); jj = ceil((n-1)/ii); kk = 1;
            figure(1); clf;
            for i = 1:ii
                for j = 1:jj
                    subplot(ii,jj,kk); hold on;
                    plot(X0 (:,1),X0 (:,kk+1),'k.',MS{:});
                    plot(Xr (:,1),Xr (:,kk+1),'ko',MS{:});
                    plot(Xf (:,1),Xf (:,kk+1),'g.',MS{:});
                    plot(Ff (:,1),Ff (:,kk+1),'b*',MS{:});
                    plot(Fbf(:,1),Fbf(:,kk+1),'r*',MS{:}); box on; hold off;
                    xlabel(VNAME{1},FS{:}); ylabel(VNAME{kk+1},FS{:});
                    if kk==1
                        if rm
                            legend('included data','removed data','included fits','removed fits','EM prev. guess','EM best fit',FS{:},'Location','northoutside');
                        else
                            legend('data','fits','EM prev. guess','EM best fit',FS{:},'Location','northoutside');
                        end
                    end
                    kk = kk+1;
                    if (kk>n-1); break; end
                end
            end
            sgtitle(['it ',num2str(it),';  misfit = ',num2str(bestfit,4),';  removed ',num2str(length(Ir)),' outliers;'])
            drawnow
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
        % determine outlier                        % add tiny random number to avoid zero criterion
        minA      = min(min(0,Af+negtol),[],2)     + rand(length(Af),1).*1e-16;  % row-wise lowest neg. mix. prop.
        maxR      = max(abs(Xf-X0)./mean(X0),[],2) + rand(length(Af),1).*1e-16;  % row-wise highest data misfit
        ir        = find(minA==min(minA)).*(abs(min(minA))>=max(maxR)) ... % select lowest negative mix. prop.
                  + find(maxR==max(maxR)).*(abs(min(minA))< max(maxR));    % or highest data misfit
        ir        = ir(1);
        
        % remove outlier and update removal diagnostics
        Abf(ir,:) = [];
        As        = sort(Abf);
        Xbf       = Abf*Fbf;
        Xr        = [Xr;X0(ir,:)];
        X0(ir,:)  = [];
        Ir        = [Ir;Ii(ir)];
        Ii(ir)    = [];
        rm        = length(Ir);
        
        % stop optimisation if >10% data identified as outliers
        if rm > m/10; error('Optimisation failed: too many samples identified as outliers.'); end
        
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