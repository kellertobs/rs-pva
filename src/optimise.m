% optimise: rs-pva subroutine to optimise for bestfit k-endmember compositions 
%           and corresponding mixing proportions and fitted data.
%
%   [Abf,Fbf,Xbf,DGN] = optimise(A0,F0,X0,DGN,VNAME)
%
%   A0     : input starting guess mixing proportions scaled up to measurement units (row sum = 1)
%   F0     : input starting guess endmember compositions scaled up to measurement units (row sum = 100%)
%   X0     : input starting guess sample compositions in measurement units (row sum = 100%)
%   Abf    : output best fit mixing proportions scaled up to measurement units (row sum = 1)
%   Fbf    : output best fit endmember compositions scaled up to measurement units (row sum = 100%)
%   Xbf    : output best fit sample compositions in measurement units (row sum = 100%)
%   DGN    : input/output rs-pva diagnostics structure for various data format and statistical metrics
%
% created  : 2020-03-21  Tobias Keller, University of Glasgow
% based on : 1994-09-12  Glenn Johnson, University of Utah
% license  : GNU General Public License v3.0


function [Abf,Fbf,Xbf,DGN] = optimise(A0,F0,X0,DGN,VNAME)

% initialise variables and parameters
Xbf     = X0;
Xf      = X0;
Af      = A0;
Ff      = F0;
Ir      = DGN.Ir;
Ii      = DGN.Ii;
rm      = DGN.rm;
misfit  = 1;
bestfit = 1;
Abf     = Af;
Fbf     = Ff;
it      = 1;
rng(15);
figure(3); clf;

% adjust convergence tolerances and step parameters
dft = [1e-6, 0.01, 1e4, 0.5, 0.0];
par = input(['->  Adjust convergence parameters as list: [msftol,negtol,maxits,alpha] \n' ...
             '    msftol: misfit tolerance for convergence          (dft = ',num2str(dft(1)),') \n' ...
             '    negtol: tolerance for negative mixing proportions (dft = ',num2str(dft(2)),') \n' ...
             '    maxits: unsuccesful iterations to outlier removal (dft = ',num2str(dft(3)),') \n' ...
             '    alpha : step size for iterative update            (dft = ',num2str(dft(4)),') \n' ...
             '    beta  : step size for randomised update           (dft = ',num2str(dft(4)),') \n']);
if isempty(par); par = dft; end

msftol  = par(1);
negtol  = par(2);
maxits  = par(3);
alpha   = par(4);
beta    = par(5);

% set preference for plotting progress
dft = 1;
plt = input('->  Every how often do you wish to plot progress of optimisation (enter # it; dft=1; no plots=0)? \n');
if isempty(plt); plt = dft; end

% initialise residuals
res_Adf = (X0(Ii,:)*X0(Ii,:).')/(Ff*X0(Ii,:).') - Af(Ii,:);
res_Ann = max(-negtol,Af(Ii,:)) - Af(Ii,:);

% iterative optimisation algorithm, minimise misfit function below tolerance
while misfit > msftol && it < maxits
    
    % update mixing proportions (no smaller than tolerance, unity sum)
    Af(Ii,:) = Af(Ii,:) + alpha .* res_Adf;
    Af(Ii,:) = Af(Ii,:) + alpha .* res_Ann;
    Af(Ii,:) = Af(Ii,:) + beta  .* (normrnd(Abf(Ii,:),abs(Abf(Ii,:)).*min(1,bestfit)) - Af(Ii,:));
    Af = Af./sum(Af,2);
    
    % update EM compositions (no negatives, 100% sum)
    Ff = (Af(Ii,:).'*Af(Ii,:))\(Af(Ii,:).'*X0(Ii,:));  % Lsq-solve for EM comp
    Ff = max(1e-3,Ff);  % zero out negative EM comp
    Ff = Ff./sum(Ff,2).*100;  % normalise EM comp to 100% sum
        
    % update data fit
    Xf(Ii,:) = Af(Ii,:)*Ff;  
    
    % update residuals
    res_Adf = (X0(Ii,:)*X0(Ii,:).')/(Ff*X0(Ii,:).') - Af(Ii,:);
    res_Ann = max(-negtol,Af(Ii,:)) - Af(Ii,:);
    
    % update residual norms
    msft_dtft = norm(res_Adf,2)./norm(Af(Ii,:),2);  % calculate data misfit 
    msft_ngpr = norm(res_Ann,2)./norm(Af(Ii,:),2);  % calculate neg. mix. prop. misfit
    misfit    = msft_dtft + msft_ngpr;  % get combined misfit
    
    if ~mod(it,plt) 
        % visualise best fit k-EM model and fitted data
        DGN.Ii = Ii; DGN.Ir = Ir; DGN.rm = rm;
        visualise({X0,Xf,Fbf},{'data','fitted data','fitted EM'},['it ',num2str(it),';  misfit = ',num2str(misfit,4),';  removed ',num2str(length(Ir)),' outliers;'],DGN,VNAME)
        figure(3); 
        plot(it,log10(msft_dtft),'bo',it,log10(msft_ngpr),'gd',it,log10(misfit),'r^','MarkerSize',5,'LineWidth',1.5); hold on; box on;
        if it==plt
            title('Convergence of Optimisation','FontSize',16); 
            xlabel('iterations','FontSize',14); ylabel('log10 misfit','FontSize',14);
            set(gca,'LineWidth',1.5);
        end
        
        % report new best fit and diagnostics
        if                  it <  10
            disp(['    ---  it    ',int2str(it),';  misfit ',num2str(misfit,'%1.3e'),';  removed ',int2str(rm),' outliers;'])
        elseif it >= 10  && it < 100
            disp(['    ---  it   ',int2str(it),';  misfit ',num2str(misfit,'%1.3e'),';   removed ',int2str(rm),' outliers;'])
        elseif it >= 100 && it < 1000
            disp(['    ---  it  ',int2str(it),';  misfit ',num2str(misfit,'%1.3e'),';  removed ',int2str(rm),' outliers;'])
        else
            disp(['    ---  it ',int2str(it),';  misfit ',num2str(misfit,'%1.3e'),';  removed ',int2str(rm),' outliers;'])
        end
    end

    % update best fit 
    if misfit < bestfit || it == 1
        % update best fit solution
        bestfit = misfit;
        Fbf     = Ff;
        Abf     = Af;
        Xbf     = Xf;
    end
    if misfit > 1e3*bestfit
        error(['  ! Optimisation algorithm diverged after ',int2str(it),' iterations.']);
    end
   
    it = it+1; % increment iteration count
    if it == maxits
        disp(['  ! Optimisation stopped after ',int2str(maxits),' iterations; final misfit ',num2str(misfit)]);
    end
end

% bundle diagnostics for output
DGN.it = it;
DGN.bf = bestfit;
DGN.rm = rm;
DGN.Ir = Ir;
DGN.Ii = Ii;

end