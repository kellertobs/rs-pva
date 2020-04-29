
% rs_pva:  Randomised Search Polytopic Vector Analysis for Geochemical Data 
%          based on an algorithm provided by Glenn Johnson, University of Utah,  
%
%          Johnson, G.W., Ehrlich, R., Full, W., Ramos, S., 2015. 
%          Principal Components Analysis and Receptor Models in Environmental 
%          Forensics. In: Murphy, B. L., Morrison, R.D. (Eds.), Introduction 
%          to Environmental Forensics, pp. 609?653. ISBN: 9780124046962.
%
% This script performs principle component and polytopic vector analysis to
% resolve the number and composition of endmembers to resolve a mixing model 
% for a set of major and trace element geochemical analyses where mixing 
% relationships are hypothesised for a genetic model.
%
% In the first step, the script analyses the descriptive statistics and
% principal components of the data set to provide metrics by which to select
% an optimal number of endmembers and reject outliers from the data set.
%
% In the second step, the script optimises the endmember compositions and
% mixing proportions to fit each sample in the data set. The routine in this
% script combines some of the steps of the original algorithm with a randomised 
% search algorithm to determine the best fit endmember compositions and mixing 
% proportions and reject further outlier points. The solution is found by
% minimising a misfit function that penalises negative mixing proportions
% and optimises data fit while requiring endmember and fitted sample
% compositions to have positive values and unit sum.
% 
% The script takes as input tabulated data placed in the sub-directory 'data'. 
% Input data should be formatted as text files (*.csv, *.txt, *.dat), with 
% analyses given in [wt %] with one row per sample, one column per variable 
% (species, oxide, etc.), with project tag or dataset name in first entry on 
% top left, variable names in first row, and sample names/numbers in first 
% column, for example:
%
% EXAMPLE DATA TABLE:
% 
% Example_Data  VAR_1   VAR_2   VAR_3   ...
% SMP_01         1.00    2.00    3.00   ...
% SMP_02         5.00    4.00    3.00   ...
% SMP_03        10.00   20.00   30.00   ...
% ...            ...     ...     ...    ...
%
% Alternatively, the script takes input files in Matlab format (.mat) with
% data places in array 'X' of dimensions (# samples x # variables), project 
% tag or dataset name in string 'PRJCT', and sample and variable names in 
% cell arrays SNAME, and VNAME, respectively. If data is generated for
% testing purposes, use array 'Xt' to provide true data values, 'At' for
% true mixing proportions, and 'Ft' for true endmember compositions.
%
% created  : 2020-03-21  Tobias Keller, University of Glasgow
% based on : 1994-09-12  Glenn Johnson, University of Utah
% license  : GNU General Public License v3.0


clear variables;       % clear workspace
warning('off','all');  % turn off warning
addpath ./src ./data   % add paths to source code and data directories

%*****  PREP: LOAD DATA  **************************************************

% set data file to read
dft  = 'DATA';
file = input(['->  What dataset would you like to process? \n', ...
              '    Give full name of datafile (dft = DATA): \n'],'s');
if isempty(file); file = dft; end

% read file to load dataset
if strcmp(file(end-3:end),'.csv') || ...  % load data from text file
   strcmp(file(end-3:end),'.txt') || ...
   strcmp(file(end-3:end),'.dat')
    DATA  = readtable(file);
    PRJCT = DATA.Properties.VariableNames(1    ); PRJCT = PRJCT{:};
    VNAME = DATA.Properties.VariableNames(2:end);
    SNAME = DATA{:    ,1};
    X     = DATA{:,2:end};
else % assume data is provided as Matlab file
    load(file);
end

% store number of samples and variables in diagnostics structure
[m,n] = size(X);
DGN.m = m; DGN.n = n;
DGN.Ii     = (1:m).';
DGN.Ir     = [];

% plot unprocessed data
if exist('Ft','var')  % if true EM known
    visualise({X,Xt,Ft},{'data','true comp.','true EM'},'Unprocessed Data',DGN,VNAME)
else
    visualise({X},{'data'},'Unprocessed Data',DGN,VNAME)
end

%*****  STEP 1: ANALYSE DATA  *********************************************

disp(' ');
disp(['Analyse dataset: ',PRJCT,' (',int2str(m),' samples, ',int2str(n),' variables):']);
disp(' ');

% analyse dataset
[Xns,Xnr,Xsq,F,A,Fpp,App,DGN] = analyse(X,DGN);

% set tolerances for EM selections
dft = [0.6,98,0.05];
tol = input(['->  Adjust endmember selection tolerances as list [CDtol,CVtol] \n' ...
             '    CDtol: selection tolerance for coefficients of data fit (dft = ',num2str(dft(1)),') \n' ...
             '    CVtol: selection tolerance for cumulative variance      (dft = ',num2str(dft(2)),') \n' ...
             '    ORtol: selection tolerance for outlier data removal     (dft = ',num2str(dft(3)),') \n' ]);
if isempty(tol); tol = dft; end
CDtol  = tol(1);
CVtol  = tol(2);
ORtol  = tol(3);

% select number of EMs 'k'
kCD    = 1 + find(min(DGN.CDvar  )>CDtol,1);
kCV    =     find(    DGN.SVD(:,4)>CVtol,1);
k      = max(kCD,kCV);
DGN.k  = k;

% do varimax rotation and initial outlier removal for chosen EMs
repeat = 1;
while repeat
    % truncate transformed EM comp 'Fpp', mix prop 'App', scaled data array 'Q'
    Appk = App(:,1:k);
    Fppk = Fpp(1:k,:);
    Qsqk = Appk*Fppk;
    
    % do varimax of mix prop and calculate EM comp
    Appvm = varimax(Appk);  % Returns scaled varimax mix prop Appvm
    
    % calculate scaled EM comp Fppvm
    Fppvm = (Appvm'*Appvm)\(Appvm'*Qsqk);

    % scale back EM comp 'Fvm' and mix prop 'Avm' to measurement units
    [Avm,Fvm] = scaleup(Fppvm,Appvm,Xns);
    
    % get fitted, row-sum normalised data 'X0' in for k-EM model
    X0 = Avm*Fvm;
    
    % remove outlier points with bad fit to k-EM model
    ir         = find(std(X0-Xns,[],2).^2./std(Xns(:)).^2 > ORtol);
    DGN.Ir     = ir;
    DGN.Ii     = (1:m).';
    DGN.Ii(ir) = [];
    DGN.rm     = length(DGN.Ir);
    
    % visualise principle components and fitted data for k-EM model
    PC = eye(k);
    PCNAME = {}; for i=1:k; PCNAME = [PCNAME,{['PC',int2str(i)]}]; end
    visualise({Appvm,PC},{'transf. data','principal comp.'},['Transformed data in ',num2str(k),' principle component space;'],DGN,PCNAME)
    
    % plot selection metrics
    FS = {'FontSize',14}; MS = {'MarkerSize',8}; LW = {'LineWidth',1.5};
    figure(2); clf;
    sgtitle(['Selected ',int2str(k),' endmembers for mixing model'],FS{:})
    subplot(1,3,1)
    plot(DGN.SVD(:,1),DGN.SVD(:,2),'k-',DGN.SVD(:,1),DGN.SVD(:,2),'ko',MS{:},LW{:}); hold on; box on; axis tight;
    plot(DGN.SVD(k,1),DGN.SVD(k,2),'ro',MS{:},LW{:});
    set(gca,LW{:});
    title('log10 singular values',FS{:});
    subplot(1,3,2)
    plot([DGN.SVD(1,1),DGN.SVD(end,1)],[CVtol,CVtol],'b-',MS{:},LW{:}); hold on; box on; axis tight;
    plot(DGN.SVD(:,1),DGN.SVD(:,4),'k-',DGN.SVD(:,1),DGN.SVD(:,4),'ko',MS{:},LW{:});
    plot(DGN.SVD(k,1),DGN.SVD(k,4),'ro',MS{:},LW{:});
    set(gca,LW{:});
    title('cum. variance',FS{:})
    subplot(1,3,3)
    plot([DGN.SVD(1,1),DGN.SVD(end,1)],[CDtol,CDtol],'b-',MS{:},LW{:}); hold on; box on; axis tight;
    plot(DGN.SVD(:,1),[nan,min(DGN.CDvar)].','k-',DGN.SVD(:,1),[nan,min(DGN.CDvar)].','ko',MS{:},LW{:});
    plot(DGN.SVD(k,1),min(DGN.CDvar(:,k-1)),'ro',MS{:},LW{:});
    set(gca,LW{:});
    title('min. CD Coefficients',FS{:})
    
    % decide to procede with selected model or change selection
    dft  = 1;
    prcd = input(['->  Proceed with ',num2str(k),'-EM model (1, dft), adjust number of EM (>1) or start over (0)? \n']);
    if isempty(prcd); prcd = dft; end
    
    if ~prcd;    return;     end
    if  prcd>1;  k = prcd;   end
    if  prcd==1; repeat = 0; end
    DGN.k = k;
end


%*****  STEP 2: PERFORM PVA  **********************************************

% Choose method to initialise polytope
dft         = 1;
init_method = input(['->  Choose PVA initialisation method (dft = 1): \n' ...
                     '    1 - EXTREME   (Initialise on k mutually extreme samples) \n' ...
                     '    2 - FUZZY     (Initialise on k fcm centroids, requires Matlab Fuzzy Toolbox) \n' ...
                     '    3 - PRINCIPAL (Initialise on k varimax principal components) \n']);
if isempty(init_method); init_method = dft; end

if     init_method == 1  % Initialise on k mutually extreme samples (Full, et al., 1981)
    [App0,Fpp0] = extreme(Appvm,Fppvm,Xns,DGN);
    [~,F0]      = scaleup(Fpp0,App0,Xns);
    F0 = max(0,F0); F0 = F0./sum(F0,2).*100;
    A0 = (X0*X0.')/(F0*X0.');
elseif init_method == 2  % Initialise on k fcm centroids (Full, et al., 1982)
    disp(' ');
    [F0,~,~] = fcm(X0,k);
    F0 = max(0,F0); F0 = F0./sum(F0,2).*100;
    A0 = (X0*X0.')/(F0*X0.');
    disp(' ');
elseif init_method == 3  % Initialise on k varimax principal components
    F0 = Fvm;
    F0 = max(0,F0); F0 = F0./sum(F0,2).*100;
    A0 = (X0*X0.')/(F0*X0.');
end

% visualise initial polytope and fitted data for k-EM model
visualise({Xns,X0,F0},{'data','fitted data','EM initial guess'},['Fitted data with ',num2str(k),'-EM starting guess;',int2str(DGN.rm),' outliers removed.'],DGN,VNAME)

% optimise for best fit EM compositions 'Fbf', mixing proportions 'Abf', and fitted data 'Xbf'
[Abf,Fbf,Xbf,DGN] = optimise(A0,F0,X0,DGN,VNAME);

% report final data fit and EM compositions
DGN.CD = (std(Xns).^2-std(Xbf-Xns).^2)./std(Xns).^2;
[~,i]  = sort(Fbf(:,1));
disp(' ');
disp('==============================================================');
disp(' ');
disp(['==> CONVERGED: ',int2str(DGN.it),' iterations; final misfit ',num2str(DGN.bf,3),'; ',int2str(DGN.rm),' outliers removed.']);
disp(' ');
fprintf(1,'                         '); 
for nn=1:n
    fprintf(1,'%s   ',VNAME{nn});
    for ll = 1:5-length(VNAME{nn}); fprintf(1,' '); end
end; fprintf(1,'\n');
disp([    '==> FINAL FIT: VAR. CD = ',num2str(DGN.CD                 ,'%1.3f   ')]);  
disp([    '               NORM CD = ',num2str(norm(DGN.CD,2)./sqrt(n),'%1.3f   ')]);    fprintf(1,'\n');
disp(     '==> FINAL EM COMPOSITIONS:');
disp(' ');
for nn=1:n
    fprintf('%s  ',VNAME{nn}); 
    for ll = 1:8-length(VNAME{nn}); fprintf(1,' '); end
    for kk=1:k; if Fbf(i(kk),nn)<10; fprintf(1,' '); end; fprintf(1,'%2.3f  ',Fbf(i(kk),nn)); end
    fprintf(1,'\n'); 
end
disp(' ');

% visualise resolved EM components and fitted data
if exist('Ft','var')
    visualise({Xns,Xbf,Fbf,Ft},{'data','best fit data','best fit EM','true EM'},['Best fit ',num2str(k),'-EM model: CD norm = ',num2str(norm(DGN.CD,2)./sqrt(n),3),', ',int2str(DGN.rm),' outliers removed'],DGN,VNAME)
else
    visualise({Xns,Xbf,Fbf},{'data','best fit data','best fit EM'},['Best fit ',num2str(k),'-EM model: CD norm = ',num2str(norm(DGN.CD,2)./sqrt(n),3),', ',int2str(DGN.rm),' outliers removed'],DGN,VNAME)
end

% end of script