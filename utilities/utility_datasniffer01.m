% A collection of simple jobs operating on results as produced by the 
% tslbatch routines. They are intended to 
% i. accomplish simple things like plotting dose response graphs
% ii. sniff out dependencies of drug effects on base line properties 

clear

% ******************* START USER INPUT ************************************
% ******************* START USER INPUT ************************************

% mat file name (without extension)
mFn='d:\hh\projects\ctx_propoSevo\data_extra\PropoSevo_ctx_wtb6_SPX_sum_raw.mat';
curIpVal=[0 .2 .5 .7];
% mat file name (without extension)
% mFn='d:\hh\projects\STN\resultsnfigs\stn_mdma.mat';
% mFn='d:\hh\projects\STN\resultsnfigs\stn_mdma_spiel.mat';
% mFn='D:\hh\projects\STN\resultsnfigs\stn_mdma_mdl_250_123.mat';
% mFn='d:\hh\projects\STN\resultsnfigs\stn_mdma_cit_new.mat';
% mFn='d:\hh\projects\STN\resultsnfigs\stn_mdma_desi_new.mat';
% mFn='d:\hh\projects\STN\resultsnfigs\stn_mdma_nom_new.mat';
% curIpVal=[0 2.5 5 10 20];

% mFn='d:\_data\otc_ctx\NO711\tslbatchPlots\NO711\tslbatch_NO711.mat';
% curIpVal=[0 10 50 250 500 1000];

% mFn='d:\_data\otc_ctx\SNAP5114\tslbatchPlots\tslbatch_SNAP.mat';
% curIpVal=[0 2.5 5 10 20];

% mFn='d:\_data\otc_ctx\muscimol\tslbatchplots\tslbatch_muscimol.mat';
% curIpVal=[0 .125 .25 .5];
% curIpVal=[0 .25];

% mFn='D:\_data\otc_ctx\AChDiazepam\Experiment1\resultsnfigs_spx\AchDiazepam_SPX.mat';
% mFn='d:\_data\otc_ctx\AChDiazepam\Experiment1\resultsnfigs_fp\AchDiazepam_FP.mat';

% % hh versions 
% mFn='d:\_data\otc_ctx\ACh\BlockDia\tslbatchPlots_BlockDia\tslbatch_BlockDia.mat';
% mFn='d:\_data\otc_ctx\ACh\AChDia\tslbatchPlots_AChDia\tslbatch_AChDia.mat';
% 
% % LL versions
% mFn='d:\_data\otc_ctx\ACh\BlockZolp\tslbatchPlots_BlockZolp\Blocker_Zolpi.mat';
% mFn='d:\_data\otc_ctx\ACh\AChZolp\tslbatchPlots_AChZolp\ACh_Zolpi.mat';


% mFn='d:\_data\otc_ctx\_oxygen\tslbatchPlots\tslbatch_oxygen.mat';
% curIpVal=1:5;
% 
% mFn='e:\_data\otc_ctx\diazepam\tslbatchPlots_dia\tslbatch_dia.mat';
% curIpVal=[0 250 1000];

% levels of independent variable (=concentrations in most cases) to deal
% with IN ASCENDING ORDER

% curIpVal=[0 50 250];
% curIpVal=[0 1 2 3]-1;

% curIpVal=[0];


% curIpVal=[0:2 4];
% curIpVal=[0 5];
% curIpVal=[0 1 2 3 4];
% curIpVal=[0 1];
% the independent var value to which data will be normalized (set to [] for
% no normalization)
normIpVal=0;
% normIpVal=[];
% normIpVal=1;


% average electrodes from cultures ??§§!!???++*##'""??
doAvEl=true;

% outliers can be eliminated in either of two ways: 
% - if outlierType is set to 'abs' outlierThresh MUST contain TWO values in
% ASCENDING order, namely the lower and upper limit of acceptable values
% - if outlierType is set to 'std' outlierThresh MUST contain ONE value,
% the logic being that for each independent value (concentration) data
% points beyond [mean +/- outlierThresh*standard deviation of the whole population]
% will be kicked
% In case you don't want to eliminate any data points set outlierType to
% 'none'
outlierType='std';
outlierThresh=2.5;

% outlierType='abs';
% outlierThresh=[0 3];
% 
outlierType='none';

% parameter(s) of interest

par={'mdSilentPerLen'};
% par={'maxSilentPerLen'};
% par={'mdBurstLen'};
% par={'asBurstLen'};
% par={'maxBurstLen'};

% par={'burstRate'};
% par={'burstIntegral'};

% par={'specpower'};
% par={'mnBurstStrength_norm'};
% par={'mnBurstPeak'};
% par={'relTimeInBurst'};

% par={'eventRate_norm'};
par={'eventRate'}; 
% par={'cvISI'}; 
% par={'eventRateInBurst'};
% par={'fractionEvOutBurst'};
% par={'fractionEvInBurst'};
% par={'mnPethEarlyFr'};
% par={'mnPethIntermedFr'};
% par={'mnPethLateFr'};

% par={'mdSTA'};


printas=[];  % don't save graphics
% printas='-djpeg97'; % save as jpeg ('working style')
% printas='-dpsc2'; % save as color post script (nice finish)

% **** job to be performed & specific options ***
job='fit';
job='dose response graph';
% job='effect dependence';

switch job
  case 'fit';
    jobOpt={'Hill'};
    jobOpt={'Hillbilly'};
    jobOpt={'linear'};

  case 'dose response graph'
    % create an unfitted and/or 2D dose response graph
    % By default, unmatched samples will be plotted in a line graph. If a
    % '.._norm' parameter is chosen the data are by definition matched. For
    % options, use any combination of the options
    %       {'matched samples','bar','list'};
    % ** note: if two parameters are specified a 2D plot will be created
    % this will kick any recording in which there is no data for any of the
    % independent values
    jobOpt={'matched samples'};
    jobOpt={''};
    % jobOpt={'list'};
    
  case 'effect dependence'
    % the parameter of interest (variable 'par' defined above) will be
    % plotted versus some property of the recording. The 'property' may be
    % any of the precomputed dependent parameters or something else (see
    % below) at a (drug) condition of your choice.
    % Choose either a precomputed dependent parameter like 'mnBurstLen' or
    % 'burstiness'
    jobOpt={'mdBurstLen'};
    % the (drug) condition of the property in question
    curAbscissIpVal=0;
end

% graphics options
labelscale('fontSz',12,'scaleFac',.8,'lineW',2,'markSz',9);

ornt='portrait';
% name & full path to figure to be saved
% figName=['d:\hh\projects\ctx_gabaUptake\rawFig\' par{1} '_' job];
figName=['d:\hh\projects\STN\paper\' par{1} '_' job];


% ******************* END USER INPUT **************************************
% ******************* END USER INPUT **************************************
% ******************* END USER INPUT **************************************
% ******************* END USER INPUT **************************************

etslconst;
nIpVal=numel(curIpVal);
% load *.mat data file
load(mFn);

% in most cases, the independent parameter is the concentration of the
% substance in question. This command plots its levels of on screen
disp('available levels of indep par:');
indepParLevel'
disp('chosen levels of indep par:');
curIpVal

% *** columnIx is the index to columns of interest in all variables
[nix,columnIx]=intersect(indepParLevel,curIpVal(:));
if ~isequal(nix,curIpVal(:))
  errordlg({'variable ''curIpVal'' contains values that do not exist in the data.',...
    'Available values: ', num2str(indepParLevel)})
  return
end

% a tribute to changes of variables across versions of tslbatch
if ~exist('expChanName','var')
  expChanName=expChanID;
end

if strcmp(job,'effect dependence')
  if numel(par)>1
    error('in job ''effect dependence'' only one dependent parameter may be specified in variable ''par''; the other must be specified in variable ''jobOpt''');
  end
  % if effect dependence is the job at hand make jobOpt first parameter
  % (=abscissa value)
  par=cat(2,jobOpt,par);
end

% reshape data such that for each parameter it's 2D
% - dependent pars
d=[];
for g=1:numel(par)
  tmpd=eval(par{g});
  tmpd=permute(tmpd,[3 2 1]);
  d=cat(3,d,tmpd);
  d(~isfinite(d))=nan;
end

% keep only columns of interest - now d and curIpVal match
d=d(:,columnIx,:);

% average electrodes from same culture?
if doAvEl
  figure; hold on; orient(ornt);
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
        plot(d(ix(1),goodIx),d(ix(2),goodIx),'o-');
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
  title('any key to continue');
  pause
  
end


% normalize, if requested
if ~isempty(normIpVal) && isfinite(normIpVal)
  % the column in d to which all others will be normalized
  col=find(curIpVal==normIpVal);
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
switch outlierType
  case 'abs'
    if numel(outlierThresh)~=2
      error('parameter ''outlierThresh'' must contain TWO values because ''outlierType'' is set to ''abs''');
    else
      if diff(outlierThresh)<=0
        error('parameter ''outlierThresh'' must contain two ASCENDING values');
      end
      badIx=d<outlierThresh(1) | d>outlierThresh(2);
      disp(['eliminating ' int2str(numel(find(badIx))) ' outliers']);
      d(badIx)=nan;
    end
  case 'std'
    badIx=d<repmat(d_av-outlierThresh*sd,size(d,1),1) | ...
      d>repmat(d_av+outlierThresh*sd,size(d,1),1);
    disp(['eliminating ' int2str(numel(find(badIx))) ' outliers']);
    d(badIx)=nan;
  case 'none'
    disp('elimination of outliers not requested');
  otherwise
    error('bad choice for parameter ''outlierType''');
end

% by default, all rows with at least one non-nan value are good
goodIx=any(isfinite(sum(d,3)),2);

% % a slightly more demanding criterion: there must be at least two non-nan
% % values per row
% goodIx=sum(isfinite(sum(d,3)),2)>1;

if ~isempty(strmatch('matched samples',jobOpt))
  % only rows without any nan are good
  goodIx=all(isfinite(sum(d,3)),2);
end

% delete all bad rows 
d=d(goodIx,:,:);
expChanName=expChanName(goodIx);

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
switch job
  case 'fit'
    curIpVal=repmat(curIpVal,size(d,1),1);
    curIpVal=curIpVal(:);
    d=d(:);
    goodIx=(isfinite(d));
    ds.crf=[curIpVal(goodIx),d(goodIx)];
    isNormalized= ~isempty(strfind(par{1},'_norm')) || ~isempty(normIpVal);
    if strncmp('Hill',jobOpt,4)


      % THIS IS THE DEFAULT CASE
      [ds,stats]=fitmaster(ds,jobOpt{1},'errType','sd','isNormalizedData',isNormalized,'invertVal',1);
%       %  THIS IS FOR STN DATA
%       [ds,stats]=fitmaster(ds,jobOpt{1},'errType','sd','isNormalizedData',isNormalized,'invertVal',nan);
    else
      [ds,stats]=fitmaster(ds,jobOpt{1},'errType','ci','isNormalizedData',isNormalized);
    end
    yl=get(gca,'ylim');
    set(gca,'ylim',[0 yl(2)]);
    orient(ornt);
    
  case 'dose response graph'
    if ~isempty(strmatch('list',jobOpt))
      clc
      [strvcat(expChanName) repmat('   ',[size(d,1) 1])  num2str(d)]
    end
    % plot
    figure; hold on; orient(ornt);
    if numel(par)<=1
      if ~isempty(strmatch('bar',jobOpt))
        bh=bar(curIpVal,d_av,'k');
        ebh=errorbar(curIpVal,d_av,d_var,'k');
        set(ebh,'marker','none','linestyle','none');
      else
        % simple plot
        if ~isempty(strmatch('matched samples',jobOpt))
          ph=plot(curIpVal,d,'ko-');
        else
          ph=plot(curIpVal,d,'ko');
        end          
        set(ph,'color',[.7 .7 .7],'linewidth',1);
        ebh=errorbar(curIpVal,d_av,d_var,'ko-');
        set(ebh,'markerfacecolor','k');
      end
      set(gca,'xtick',curIpVal);
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
      xlabel(par{1})
      ylabel(par{2});
      legend(ph,num2str(indepParLevel(columnIx)));
      nicexyax
    end
    
  case 'effect dependence'
    % the independent value of the abscissa par 
    abscissParCol=find(curIpVal==curAbscissIpVal);
%     % if y values to be plotted include values at the independent par
%     % value used for normalization, kick them, because their y values will
%     % all be 1 anyways
%     killColIx=columnIx==indepParNormIx;
%     d(:,killColIx,:)=nan;
    xlab=[jobOpt{1} ' at ' num2str(curIpVal(abscissParCol))];
    figure(1); clf; hold on; orient(ornt);
    ph=plot(d(:,abscissParCol,1),d(:,:,2),'o');
    nicexyax
    xlabel(xlab)
    ylabel(par{2})
    legend(ph,num2str(indepParLevel(setdiff(columnIx,indepParNormIx))));

  case 'effect dependence_old'
    par=par{1};
    d=eval(par);
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
    abscissParCol=find(indepParLevel==curAbscissIpVal);
    switch jobOpt{1}
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
        abscissVal=eval(jobOpt{1});
        abscissVal=permute(abscissVal,[3 2 1]);
        abscissVal=abscissVal(:,abscissParCol);
        xlab=[jobOpt{1} ' at ' num2str(indepParLevel(abscissParCol))];         
    end
    figure(1); clf; hold on; orient(ornt);
    ph=plot(abscissVal,d,'o');
    nicexyax
    xlabel(xlab)
    ylabel(par)
    legend(ph,num2str(indepParLevel(setdiff(columnIx,indepParNormIx))));
  otherwise 
    error('illegal job')
end
drawnow

% ----- printing -----
if ~isempty(printas), 
  print(printas,figName,'-r400'); 
end