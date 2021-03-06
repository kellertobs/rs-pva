% visualise:  rs-pva subroutine to visualise compositional data
%
%   [] = visualise(DATA,LEGEND,TITLE,DGN,VNAME)
%
%   DATA   : input cell array containing standard-sized data (X) and/or endmember arrays (F)
%   LEGEND : input cell array containing legend strings in same sequence as data arrays in DATA
%   TITLE  : input string for global figure title
%   DGN    : input rs-pva diagnostics structure as generated by main routine
%   VNAME  : input cell array with variable names as loaded from data table
%
% created  : 2020-03-21  Tobias Keller, University of Glasgow
% license  : GNU General Public License v3.0


function    [] = visualise(DATA,LEGEND,TITLE,DGN,VNAME)

FS = {'FontSize',14}; MT = {'.','o','*'}; MS = {'MarkerSize',8};
LW = {'LineWidth',1.5}; CL = {'Color','k','b','r','m','g','c'};
kk = length(DATA);
nn = size(DATA{1},2);
ii = ceil(sqrt(nn-1));
jj = ceil((nn-1)/ii);
pp = 1;
figure(1); clf;
sgtitle(TITLE,FS{1},FS{2}+2);
for i = 1:ii
    for j = 1:jj
        subplot(ii,jj,pp); box on; hold on;
        for k = 1:kk
            if     size(DATA{k},1) == DGN.m % is data array
                plot(DATA{k}(DGN.Ii,1),DATA{k}(DGN.Ii,pp+1),MT{1},MS{:},CL{1},CL{1+k},LW{:});
                plot(DATA{k}(DGN.Ir,1),DATA{k}(DGN.Ir,pp+1),MT{2},MS{:},CL{1},CL{1+k},LW{:});
                if ~isempty(DGN.Ir); LEGEND = {LEGEND{1:2*k-1},[LEGEND{2*k-1},' removed'],LEGEND{2*k:end}}; end
            else % is EM array
                plot(DATA{k}(:,1),DATA{k}(:,pp+1),MT{3},MS{:},CL{1},CL{1+k},LW{:});
            end
        end
        if pp==1; legend(LEGEND{:},FS{:},LW{:},'Location','northoutside','box','on'); end
        xlabel(VNAME{1},FS{:}); ylabel(VNAME{pp+1},FS{:});
        set(gca,LW{:});
        pp = pp+1;
        if (pp>nn-1); break; end
    end
end
drawnow;