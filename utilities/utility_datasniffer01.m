% A collection of simple jobs operating on results as produced by the 
% tslbatch routines. They are intended to 
% i. accomplish simple things like plotting dose response graphs
% ii. sniff out dependencies of drug effects on base line properties 

clear

% ******************* START USER INPUT ************************************
% ******************* START USER INPUT ************************************

% mat file name (without extension)
% mFn='D:\hh\projects\STN\m\figures\stn_mdma_spiel.mat';
mFn='d:\hh\projects\ctx_gabaUptake\resultsnfigs_SNAP\tslbatchPlots\tslbatch_SNAP.mat';
% mFn='D:\_data\otc_ctx\AChDiazepam\Experiment1\resultsnfigs_spx\AchDiazepam_SPX.mat';
% mFn='d:\_data\otc_ctx\AChDiazepam\Experiment1\resultsnfigs_fp\AchDiazepam_FP.mat';
mFn='d:\hh\_tbs\AChpilot.mat';


% levels of independent variable (=concentrations in most cases) to deal
% with IN ASCENDING ORDER
% curIpVal=[0 250];
% curIpVal=[0 2.5 5 10 20];
curIpVal=[0 1 2];

% parameter(s) of interest
% par={'asBurstLen','mdBurstLen'};
% par={'specpower'};
% par={'asBurstLen'};
% par={'mnBurstStrength_norm'};
% par={'mnBurstLen_norm','burstRate_norm'};
par={'mnSilentPerLen'};
% par={'mnBurstLen'};
% par={'relTimeInBurst_norm'};
% par={'eventRate_norm'};
% par={'mnBurstStrength_norm'};
% par={'burstRate_norm'};

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
    % jobOpt={'linear'};
  case 'dose response graph'
    % create an unfitted and/or 2D dose response graph
    % By default, unmatched samples will be plotted in a line graph. If a
    % '.._norm' parameter is chosen the data are by definition matched. For
    % options, use any combination of the options
    %       {'matched samples','bar','list'};
    % ** note: if two parameters are specified a 2D plot will be created
    jobOpt={' '};
  case 'effect dependence'
    % the parameter of interest (first element of variable 'par') will be
    % plotted versus some property of the recording under a drug condition
    % to be specified below. The 'property' may be another of the
    % precomputed dependent parameters or something else (see below). This
    % option is more general and it allows a computation of arbitrary
    % parameters from the collected results
    % Choose either a precomputed dependent parameter like 'mnBurstLen' or
    % 'burstiness'
    % jobOpt={'eventRate'};
    jobOpt={'maxBurstLen'};
    % jobOpt={'burstiness'};
    % jobOpt={'drug1_applicRank'};
    % the concentration at which the property in question shall be quantified
    ipVal=0;
end

% graphics options
labelscale('fontSz',12,'scaleFac',.6,'lineW',1.5,'markSz',10);
ornt='portrait';
% name & full path to figure to be saved
figName=['D:\hh\projects\ctx_gabaUptake\rawFig\' par{1} '_' job];


% ******************* END USER INPUT **************************************
% ******************* END USER INPUT **************************************

etslconst;
nIpVal=numel(curIpVal);
% load *.mat data file
load(mFn);

% in most cases, the independent parameter is the concentration of the
% substance in question. This command plots its levels of on screen
indepParLevel
% index to columns in all variables
[nix,columnIx]=intersect(indepParLevel,curIpVal);
if ~isequal(nix,curIpVal)
  errordlg({'variable ''curIpVal'' contains values that do not exist in the data.',...
    'Available values: ', num2str(indepParLevel)})
  return
end
    
% reshape data such that it's 2D
% - specific independent pars
% drug1_applicRank=permute(drug1_applicRank,[3 2 1]); 
% - dependent pars

d=[];
for g=1:numel(par)
  tmpd=eval(par{g});
  tmpd=permute(tmpd,[3 2 1]);
  d=cat(3,d,tmpd);
  d(~isfinite(d))=nan;
end

switch job
  case 'fit'
    curIpVal=repmat(curIpVal,size(d,1),1);
    d=d(:,columnIx,:);
    curIpVal=curIpVal(:);
    d=d(:);
    goodIx=(isfinite(d));
    ds.crf=[curIpVal(goodIx),d(goodIx)];
    isNormalized= ~isempty(strfind(par{1},'_norm'));
    if strmatch('Hill',jobOpt)
      if ~isNormalized
        error('Hill fit works only on normalized data');
      end
      
%       % *********** a hack for importing NO-711 data
%       load d:\hh\projects\ctx_gabaUptake\rawFig\plot_wsdata_export.mat
%       % get rid of nans
%       d(isnan(d(:,2)),:)=[]
%       % convert conc to uM
%       d(:,1)=d(:,1)/1000;
%       ds.crf=d;
%       % ********************************************
      
      [ds,stats]=fitmaster(ds,jobOpt{1},'errType','ci','isNormalizedData',isNormalized,'invertVal',1);
    else
      [ds,stats]=fitmaster(ds,jobOpt{1},'errType','ci','isNormalizedData',isNormalized);
    end
    yl=get(gca,'ylim');
    set(gca,'ylim',[0 yl(2)]);
    orient(ornt);
    
  case 'dose response graph'
    d=d(:,columnIx,:);
    if ~isempty(strmatch('matched samples',jobOpt))
      % delete all rows with any nan
      d=d(all(isfinite(sum(d,3)),2),:);
    end
    % mean, sd, ci
    d_av=nanmean(d);
    sd=nanstd(d);
    finIx=isfinite(d);
    n=sum(finIx);
    sem=sd./sqrt(n);
    ci=sem.* -permute(tinv(.05/2,permute(n,[3 2 1])-1),[3 2 1]);
    d_var=ci;    
    % d_var=sem;
    if ~isempty(strmatch('list',jobOpt))
      clc
      [strvcat(expChanID) repmat('   ',[size(d,1) 1])  num2str(d)]
    end
    % plot
    figure(1); clf; hold on; orient(ornt);
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
    par=par{1};
    if isempty(strfind(par,'_norm'))
      par=[par '_norm'];
    end
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
    col=find(indepParLevel==ipVal);
    switch jobOpt{1}
      case 'burst length'
        abscissVal=permute(mdBurstLen,[3 2 1]);
        abscissVal=abscissVal(:,col);
        xlab=['median of burst length at ' num2str(indepParLevel(col))]; 
      case 'variability burst length'
        abscissVal=permute(cvBurstLen,[3 2 1]);
        abscissVal=abscissVal(:,col);
        xlab=['CV of burst length at ' num2str(indepParLevel(col))]; 
      case 'asymmetry burst length'
        abscissVal=permute(asBurstLen,[3 2 1]);
        abscissVal=abscissVal(:,col);
        xlab=['asymmetry of burst length at ' num2str(indepParLevel(col))]; 
      case 'burstiness'
        abscissVal=repmat(nan,[size(Tsl,1) 1]);
        for g=1:size(Tsl,1)
          tmpIsi=diff(Tsl{g,col});
          abscissVal(g)=std(tmpIsi)./mean(tmpIsi);
        end
        % coefficient of variation of inter-spike interval
        xlab=['CV of ISI at ' num2str(indepParLevel(col))]; 
      otherwise
        abscissVal=eval(jobOpt{1});
        abscissVal=permute(abscissVal,[3 2 1]);
        abscissVal=abscissVal(:,col);
        xlab=[jobOpt{1} ' at ' num2str(indepParLevel(col))];         
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