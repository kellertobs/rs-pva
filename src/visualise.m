function    [] = visualise(DATA,LEGEND,TITLE,DGN,VNAME)

FS = {'FontSize',14}; MT = {'.','o','*'}; MS = {'MarkerSize',8};
LW = {'LineWidth',1.5}; CL = {'Color','k','b','r','m','g','c'};
ii = ceil(sqrt(DGN.n-1));
jj = ceil((DGN.n-1)/ii);
kk = length(DATA);
pp = 1;
figure(1); clf;
sgtitle(TITLE,FS{1},FS{2}+2);
for i = 1:ii
    for j = 1:jj
        subplot(ii,jj,pp); box on; hold on;
        for k = 1:kk
            if     length(DATA{k}) == DGN.m % is data array
                plot(DATA{k}(DGN.Ii,1),DATA{k}(DGN.Ii,pp+1),MT{1},MS{:},CL{1},CL{1+k},LW{:});
                plot(DATA{k}(DGN.Ir,1),DATA{k}(DGN.Ir,pp+1),MT{2},MS{:},CL{1},CL{1+k},LW{:});
                if ~isempty(DGN.Ir); LEGEND = {LEGEND{1:2*k-1},[LEGEND{2*k-1},' removed'],LEGEND{2*k:end}}; end
            elseif length(DATA{k}) == DGN.n % is EM array
                plot(DATA{k}(:,1),DATA{k}(:,pp+1),MT{3},MS{:},CL{1},CL{1+k},LW{:});
            end
        end
        if pp==1; legend(LEGEND{:},FS{:},LW{:},'Location','northoutside','box','on'); end
        xlabel(VNAME{1},FS{:}); ylabel(VNAME{pp+1},FS{:});
        set(gca,LW{:});
        pp = pp+1;
        if (pp>DGN.n-1); break; end
    end
end
drawnow;