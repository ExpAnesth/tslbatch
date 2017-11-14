function varargout=datasniffer_tslbatch_01(ds)
% A collection of simple jobs operating on results as produced by the 
% tslbatch routines. They are intended to 
% i. accomplish simple things like plotting dose response graphs
% ii. sniff out dependencies of drug effects on base line properties 

% -------------------------------------------------------------------------
% 0. PRELUDE, GLOBAL SETTINGS & USER INPUT
% -------------------------------------------------------------------------

if ~isstruct(ds)
  error('input arg ''ds'' must be a struct');
end

ornt='portrait';

% -------------------------------------------------------------------------
% 2. PREPS & CHECKING SETTINGS
% -------------------------------------------------------------------------

% name of matfile as part of title
mFn=strfind(ds.mFnFull,'\');
if ~isempty(mFn) && mFn(end)<numel(ds.mFnFull)
  mFn=ds.mFnFull(mFn(end)+1:end);
else
  mFn=ds.mFnFull;
end
% full path & name of figure file to be saved
figFn=[ds.figPath mFn '_' ds.par{1,1}];

nIpVal=numel(ds.curIpVal);

% -------------------------------------------------------------------------
% 3. READ & HONE DATA
% -------------------------------------------------------------------------

% load *.mat data file
load(mFn);

% in most cases, the independent parameter is the concentration of the
% substance in question. This command plots its levels of on screen
disp('available levels of indep ds.par:');
indepParLevel'
disp('chosen levels of indep ds.par:');
ds.curIpVal

% *** columnIx is the index to columns of interest in all variables
[nix,columnIx]=intersect(indepParLevel,ds.curIpVal(:));
if ~isequal(nix,ds.curIpVal(:))
  errordlg({'variable ''ds.curIpVal'' contains values that do not exist in the data.',...
    'Available values: ', num2str(indepParLevel)})
  return
end

% a tribute to changes of variables across versions of tslbatch
if ~exist('expChanName','var')
  expChanName=expChanID;
end

if strcmp(ds.job,'effect dependence')
  if numel(ds.par)>1
    error('in ds.job ''effect dependence'' only one dependent parameter may be specified in variable ''ds.par''; the other must be specified in variable ''ds.jobOpt''');
  end
  % if effect dependence is the ds.job at hand make ds.jobOpt first parameter
  % (=abscissa value)
  ds.par=cat(2,ds.jobOpt,ds.par);
end

% reshape data such that for each parameter it's 2D
% - dependent pars
d=[];
for g=1:numel(ds.par)
  tmpd=eval(ds.par{g});
  tmpd=permute(tmpd,[3 2 1]);
  d=cat(3,d,tmpd);
  d(~isfinite(d))=nan;
end

% keep only columns of interest - now d and ds.curIpVal match
d=d(:,columnIx,:);

% average electrodes from same culture?
if ds.doAvEl
  % figure(100); clf; hold on; orient(ornt);
  nECh=numel(expChanName);
  expName={};
  d2=nan(size(d));
  d3=[];
  ct=1;
  list=1:nECh;
  while ~isempty(list)
    curExp=expChanName{list(1)};
    commaIx=strfind(curExp,',');
    expName{ct}=curExp(1:commaIx-1);
    % index to entries of current exp
    ix=find(strncmp(curExp(1:commaIx-1),expChanName,commaIx-1));
    d2(ct,:,:)=nanmean(d(ix,:,:),1);
    % if there are two recording sites...
    if numel(ix)>1
      % (consider only first variable)
      goodIx=all(isfinite(d(ix,:,1)));
      if any(goodIx)
        % ... collect data such that the correlation between them can be
        % computed
        d3=cat(1,d3,d(ix,goodIx,1)');
        % plot(d(ix(1),goodIx),d(ix(2),goodIx),'o-');
      end
    end
    list=setdiff(list,ix);
    ct=ct+1;
  end
  % overwrite, killing the superfluous preallocated rows 
  d=d2(1:ct-1,:,:);
  
  % remove normalized control values from d3, if any
  d3(all(d3==1,2),:)=[];
  % finally, correlation
  tmp=corrcoef(d3);
  disp(['correlation between recording sites for (first) chosen parameter: ' num2str(tmp(1,2))]);
  
  axis tight
  axis square
else
  expName=expChanName;
end


% normalize, if requested
if ~isempty(ds.normIpVal) && isfinite(ds.normIpVal)
  % the column in d to which all others will be normalized
  col=find(ds.curIpVal==ds.normIpVal);
  d=d./repmat(d(:,col,:),[1 size(d,2) 1]);
  d(isinf(d))=nan;
end

% run of the mill-preprocessing
d_av=nanmean(d);
sd=nanstd(d);
finIx=isfinite(d);
n=sum(finIx);
sem=sd./sqrt(n);
ci=sem.* -permute(tinv(.05/2,permute(n,[3 2 1])-1),[3 2 1]);
d_var=sd;
% d_var=sem;

% kick outliers, if requested
switch ds.outlierType
  case 'abs'
    if numel(ds.outlierThresh)~=2
      error('parameter ''ds.outlierThresh'' must contain TWO values because ''ds.outlierType'' is set to ''abs''');
    else
      if diff(ds.outlierThresh)<=0
        error('parameter ''ds.outlierThresh'' must contain two ASCENDING values');
      end
      badIx=d<ds.outlierThresh(1) | d>ds.outlierThresh(2);
      disp(['eliminating ' int2str(numel(find(badIx))) ' outliers']);
      d(badIx)=nan;
    end
  case 'std'
    badIx=d<repmat(d_av-ds.outlierThresh*sd,size(d,1),1) | ...
      d>repmat(d_av+ds.outlierThresh*sd,size(d,1),1);
    disp(['eliminating ' int2str(numel(find(badIx))) ' outliers']);
    d(badIx)=nan;
  case 'none'
    disp('elimination of outliers not requested');
  otherwise
    error('bad choice for parameter ''ds.outlierType''');
end

% by default, all rows with at least one non-nan value are good
goodIx=any(isfinite(sum(d,3)),2);

% % a slightly more demanding criterion: there must be at least two non-nan
% % values per row
% goodIx=sum(isfinite(sum(d,3)),2)>1;

if ~isempty(strmatch('matched samples',ds.jobOpt))
  % only rows without any nan are good
  goodIx=all(isfinite(sum(d,3)),2);
end

% delete all bad rows 
d=d(goodIx,:,:);
expChanName=expChanName(goodIx);
expName=expName(goodIx);

% run of the mill-processing reloaded
d_av=nanmean(d);
sd=nanstd(d);
finIx=isfinite(d);
n=sum(finIx);
sem=sd./sqrt(n);
ci=sem.* -permute(tinv(.05/2,permute(n,[3 2 1])-1),[3 2 1]);
d_var=sd;
% d_var=sem;
  

% **** jobs proper ****
switch ds.job
  case 'list'
    clc
    disp(['*** ' ds.par{1} ':']);
    disp('')
    [char(expName) repmat('   ',[size(d,1) 1])  num2str(d)]
    
  case 'fit'
    ds.curIpVal=repmat(ds.curIpVal,size(d,1),1);
    ds.curIpVal=ds.curIpVal(:);
    d=d(:);
    goodIx=(isfinite(d));
    ds.crf=[ds.curIpVal(goodIx),d(goodIx)];
    isNormalized= ~isempty(strfind(ds.par{1},'_norm')) || ~isempty(ds.normIpVal);
    if strncmp('Hill',ds.jobOpt,4)


      % THIS IS THE DEFAULT CASE
      [ds,stats]=fitmaster(ds,ds.jobOpt{1},'errType','sd','isNormalizedData',isNormalized,'invertVal',1);
%       %  THIS IS FOR STN DATA
%       [ds,stats]=fitmaster(ds,ds.jobOpt{1},'errType','sd','isNormalizedData',isNormalized,'invertVal',nan);
    else
      [ds,stats]=fitmaster(ds,ds.jobOpt{1},'errType','ci','isNormalizedData',isNormalized);
    end
    yl=get(gca,'ylim');
    set(gca,'ylim',[0 yl(2)]);
    orient(ornt);
    
  case 'dose response graph'
    if ~isempty(strmatch('list',ds.jobOpt))
      clc
      [strvcat(expChanName) repmat('   ',[size(d,1) 1])  num2str(d)]
    end
    % plot
    figure; hold on; orient(ornt);
    if numel(ds.par)<=1
      if ~isempty(strmatch('bar',ds.jobOpt))
        bh=bar(ds.curIpVal,d_av,'k');
        ebh=errorbar(ds.curIpVal,d_av,d_var,'k');
        set(ebh,'marker','none','linestyle','none');
      else
        % simple plot
        if ~isempty(strmatch('matched samples',ds.jobOpt))
          ph=plot(ds.curIpVal,d,'ko-');
        else
          ph=plot(ds.curIpVal,d,'ko');
        end          
        set(ph,'color',[.7 .7 .7],'linewidth',1);
        ebh=errorbar(ds.curIpVal,d_av,d_var,'ko-');
        set(ebh,'markerfacecolor','k');
      end
      set(gca,'xtick',ds.curIpVal);
      nicexyax;
    else
      cm=flipud(coma('blues','ncols',nIpVal));
      set(gca,'colororder',cm);
      ph=plot(d(:,:,1)',d(:,:,2)','k-');
      set(ph,'color',[.5 .5 .5]);      
      ph=[];
      for g=1:nIpVal
        ph(g)=plot(d(:,g,1),d(:,g,2),'o');
        set(ph(g),'color',[.5 .5 .5],'markerfacecolor',cm(g,:));
      end
      xlabel(ds.par{1})
      ylabel(ds.par{2});
      legend(ph,num2str(indepParLevel(columnIx)));
      nicexyax
    end
    
  case 'effect dependence'
    % the independent value of the abscissa ds.par 
    abscissParCol=find(ds.curIpVal==ds.curAbscissIpVal);
%     % if y values to be plotted include values at the independent ds.par
%     % value used for normalization, kick them, because their y values will
%     % all be 1 anyways
%     killColIx=columnIx==indepParNormIx;
%     d(:,killColIx,:)=nan;
    xlab=[ds.jobOpt{1} ' at ' num2str(ds.curIpVal(abscissParCol))];
    figure(1); clf; hold on; orient(ornt);
    ph=plot(d(:,abscissParCol,1),d(:,:,2),'o');
    nicexyax
    xlabel(xlab)
    ylabel(ds.par{2})
    legend(ph,num2str(indepParLevel(setdiff(columnIx,indepParNormIx))));

  case 'effect dependence_old'
    ds.par=ds.par{1};
    d=eval(ds.par);
    d=permute(d,[3 2 1]);
    if exist('Tsl','var')
      Tsl=permute(Tsl,[3 2 1]);
    end
    if exist('Etsl','var')
      Etsl=permute(Etsl,[3 2 1]);
    end
    % effect: compute for all but for control values because these will be
    % 1 anyways
    tmpColIx=setdiff(columnIx,indepParNormIx);
    d=d(:,tmpColIx);
    % the column in d in which to quantify property
    abscissParCol=find(indepParLevel==ds.curAbscissIpVal);
    switch ds.jobOpt{1}
      case 'burst length'
        abscissVal=permute(mdBurstLen,[3 2 1]);
        abscissVal=abscissVal(:,abscissParCol);
        xlab=['median of burst length at ' num2str(indepParLevel(abscissParCol))]; 
      case 'variability burst length'
        abscissVal=permute(cvBurstLen,[3 2 1]);
        abscissVal=abscissVal(:,abscissParCol);
        xlab=['CV of burst length at ' num2str(indepParLevel(abscissParCol))]; 
      case 'asymmetry burst length'
        abscissVal=permute(asBurstLen,[3 2 1]);
        abscissVal=abscissVal(:,abscissParCol);
        xlab=['asymmetry of burst length at ' num2str(indepParLevel(abscissParCol))]; 
      case 'burstiness'
        abscissVal=repmat(nan,[size(Tsl,1) 1]);
        for g=1:size(Tsl,1)
          tmpIsi=diff(Tsl{g,abscissParCol});
          abscissVal(g)=std(tmpIsi)./mean(tmpIsi);
        end
        % coefficient of variation of inter-spike interval
        xlab=['CV of ISI at ' num2str(indepParLevel(abscissParCol))]; 
      otherwise
        abscissVal=eval(ds.jobOpt{1});
        abscissVal=permute(abscissVal,[3 2 1]);
        abscissVal=abscissVal(:,abscissParCol);
        xlab=[ds.jobOpt{1} ' at ' num2str(indepParLevel(abscissParCol))];         
    end
    figure(1); clf; hold on; orient(ornt);
    ph=plot(abscissVal,d,'o');
    nicexyax
    xlabel(xlab)
    ylabel(ds.par)
    legend(ph,num2str(indepParLevel(setdiff(columnIx,indepParNormIx))));
  otherwise 
    error('illegal ds.job')
end

varargout{1}=d;
varargout{2}=expName;
varargout{3}=indepParLevel;

drawnow

% ----- printing -----
if ~isempty(ds.printas), 
  print(ds.printas,'-r300',figFn); 
end


