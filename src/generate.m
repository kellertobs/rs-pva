% generate samples for testing PVA
m = 50;   % number of samples
n = 5;    % number of variables
k = 3;    % number of endmembers
a = 1;    % amplitude of added noise [wt %]

% generate true endmember compositions
Ft = rand(k,n); 
Ft = Ft./sum(Ft,2).*100;  % in wt %

% generate true mixing proportions
At = rand(m,k);
At = At./sum(At,2);  % in wt fract.

% generate true sample compositions
Xt = At*Ft;

% perturb sample compositions with normally distributed noise
X = Xt + a.*(rand(size(Xt))-0.5);

% add 5% outlier points
o       = floor(0.05*m);
io      = randi(m,o,1);
X(io,:) = X(io,:) .* (1-(randi(2,o,n)-1).*0.5);

% set sample, variable, and project names
SNAME = {}; for i = 1:m; SNAME = [SNAME;{['S',num2str(i)]}]; end
VNAME = {}; for i = 1:n; VNAME = [VNAME;{['V',num2str(i)]}]; end
PRJCT = ['Data generated ',datestr(now)]; 
PRJCT(PRJCT==' ' | PRJCT=='-') = '_';

ii = ceil(sqrt(n-1)); jj = ceil((n-1)/ii); kk = 1;
FS = {'FontSize',14}; MS = {'MarkerSize',8};
figure(1); clf;
for i = 1:ii
    for j = 1:jj
        subplot(ii,jj,kk); box on; hold on;
        plot(X (:,1),X (:,kk+1),'k.',MS{:}); 
        plot(Xt(:,1),Xt(:,kk+1),'b.',MS{:});
        plot(Ft(:,1),Ft(:,kk+1),'r*',MS{:});
        xlabel(VNAME{1},FS{:}); ylabel(VNAME{kk+1},FS{:});
        if kk==1; legend('generated data','true data','endmembers',FS{:},'Location','best'); end
        kk = kk+1;
        if (kk>n-1); break; end
    end
end
sgtitle(['Generated ',int2str(m),' samples in ',int2str(n),' variables between ',int2str(k),' endmembers'],FS{:})
        
save('../data/DATA.mat','X','Xt','At','Ft','SNAME','VNAME','PRJCT');