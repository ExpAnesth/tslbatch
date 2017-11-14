% Principal components analyis operating on results as produced by the 
% tslbatch routines. It is intended to sniff out the most dramatic changes
% between conditions

clear

% ******************* START USER INPUT ************************************
% ******************* START USER INPUT ************************************

% mat file names (without extension)
% ALL
mFn={'d:\hh\projects\ctx_ACh\_AChBlockDia\tslbatchPlots_BlockDia\Blocker_Dia_Norm1',...
  'd:\hh\projects\ctx_ACh\_AChBlockDia\tslbatchPlots_AChDia\ACh_Dia_Norm1',...
  'd:\hh\projects\ctx_ACh\_AChBlockZolp\tslbatchPlots_BlockZolp\Blocker_Zolpi_Norm1',...
  'd:\hh\projects\ctx_ACh\_AChBlockZolp\tslbatchPlots_AChZolp\ACh_Zolpi'};

% DIA ONLY
mFn={'d:\hh\projects\ctx_ACh\_AChBlockDia\tslbatchPlots_BlockDia\Blocker_Dia_Norm1',...
  'd:\hh\projects\ctx_ACh\_AChBlockDia\tslbatchPlots_AChDia\ACh_Dia_Norm1'};

% ZOL ONLY
mFn={'d:\hh\projects\ctx_ACh\_AChBlockZolp\tslbatchPlots_BlockZolp\Blocker_Zolpi_Norm1',...
  'd:\hh\projects\ctx_ACh\_AChBlockZolp\tslbatchPlots_AChZolp\ACh_Zolpi'};

% % ACh ONLY
% mFn={'d:\hh\projects\ctx_ACh\_AChBlockDia\tslbatchPlots_AChDia\ACh_Dia_Norm1',...
%   'd:\hh\projects\ctx_ACh\_AChBlockZolp\tslbatchPlots_AChZolp\ACh_Zolpi'};

% % BLOCKER ONLY
% mFn={'d:\hh\projects\ctx_ACh\_AChBlockDia\tslbatchPlots_BlockDia\Blocker_Dia_Norm1',...
%   'd:\hh\projects\ctx_ACh\_AChBlockZolp\tslbatchPlots_BlockZolp\Blocker_Zolpi_Norm1'};

% % blocker vs ACh
% mFn={'d:\hh\projects\ctx_ACh\_AChBlock\tslbatchPlots\Blocker_all',...
%   'd:\hh\projects\ctx_ACh\_AChBlock\tslbatchPlots\ACh_out'};


% levels of independent variable (=concentrations in most cases) to deal
% with IN ASCENDING ORDER
curIpVal=[1 2];
% curIpVal=[0 1];
% the independent var value to which data will be normalized (set to [] for
% no normalization)
normIpVal=1;
% normIpVal=0;
% normIpVal=[];

% the indep value at which parameters shall be analyzed
anIpVal=2;
% anIpVal=1;

% outliers can be eliminated: 
% - if outlierType is set to 'std' outlierThresh MUST contain ONE value,
% the logic being that for each independent value (concentration) data
% points beyond [mean +/- outlierThresh*standard deviation of the whole population]
% will be kicked
% In case you don't want to eliminate any data points set outlierType to
% 'none'
outlierType='std';
outlierThresh=2.5;
% outlierType='none';

% parameter(s) of interest
par={'mdBurstPeak','mdBurstStrength','burstIntegral','relTimeInBurst',...  
  'mdBurstLen','burstRate','mdSilentPerLen','specpower','eventRate',...
  'mnPethEarlyFr','mnPethIntermedFr','mnPethLateFr','mnSpPethFr'};

par={'mdBurstPeak','mdBurstStrength','burstIntegral','relTimeInBurst',...  
  'mdBurstLen','burstRate','mdSilentPerLen','eventRate',...
  'mnPethEarlyFr','mnPethIntermedFr','mnPethLateFr','mnSpPethFr'};

% LFP pars only
par={'mdBurstPeak','mdBurstStrength','burstIntegral','relTimeInBurst',...  
  'mdBurstLen','burstRate','mdSilentPerLen'};

% % SPX pars only
% par={'eventRate',...
%   'mnPethEarlyFr','mnPethIntermedFr','mnPethLateFr','mnSpPethFr'};

% number of principal components
nPC=3;

printas=[];  % don't save graphics
% printas='-djpeg97'; % save as jpeg ('working style')
% printas='-dpsc2'; % save as color post script (nice finish)

% **** job to be performed & specific options ***
job='PCA';

switch job
  case 'PCA';
    jobOpt={''};
end

% graphics options
labelscale('fontSz',10,'scaleFac',.5,'lineW',1.5,'markSz',8);

ornt='portrait';
% name & full path to figure to be saved
figName=['d:\hh\projects\ctx_ACh\rawFig\' mfilename '_' job];


% ******************* END USER INPUT **************************************
% ******************* END USER INPUT **************************************
% ******************* END USER INPUT **************************************

% =========================================================================
% ****************** PART 1: LOADING & CONCATENATING DATA ******************
% =========================================================================
nIpVal=numel(curIpVal);
nFile=numel(mFn);

[nix,anColumnIx]=intersect(curIpVal,anIpVal);
% collected (concatenated) data
collCatD=[];
% tags for the data in different files
tag=[];

for fi=1:nFile
  % load *.mat data file
  load(mFn{fi});
  
  % *** columnIx is the index to columns of interest in all variables
  [nix,columnIx]=intersect(indepParLevel,curIpVal);
  if ~isequal(nix,curIpVal)
    errordlg({'variable ''curIpVal'' contains values that do not exist in the data.',...
      'Available values: ', num2str(indepParLevel)})
    return
  end
  
  % reshape data such that for each parameter it's 2D
  d=[];
  for g=1:numel(par)
    tmpd=eval(par{g});
    % permute, picking only columns of interest
    tmpd=permute(tmpd(:,columnIx,:),[3 2 1]);
    d=cat(3,d,tmpd);
  end
  tag=cat(1,tag,repmat(fi,[size(tmpd,1) 1]));
  
  % ** d now has the following structure:
  % experiment | indep var | parameter
  
  % normalize, if requested
  if ~isempty(normIpVal) && isfinite(normIpVal)
    % the column in d to which all others will be normalized
    col=find(curIpVal==normIpVal);
    d=d./repmat(d(:,col,:),[1 size(d,2) 1]);
  end

  % get rid of infs
  d(~isfinite(d))=nan;
  
  % run of the mill-preprocessing
  d_av=nanmean(d);
  sd=nanstd(d);
  finIx=isfinite(d);
  n=sum(finIx);
  sem=sd./sqrt(n);
  ci=sem.* -permute(tinv(.05/2,permute(n,[3 2 1])-1),[3 2 1]);
  d_var=ci;
  
  % kick outliers, if requested
  switch outlierType
    case 'std'
      badIx=d<repmat(d_av-outlierThresh*sd,size(d,1),1) | ...
        d>repmat(d_av+outlierThresh*sd,size(d,1),1);
      disp(['eliminating ' int2str(numel(find(badIx))) ' outliers']);
      d(badIx)=nan;
      doReprocess=true;
    case 'none'
      disp('elimination of outliers not requested');
      doReprocess=false;
    otherwise
      error('bad choice for parameter ''outlierType''');
  end
  
  % % run of the mill-processing reloaded
  % if doReprocess
  %   d_av=nanmean(d);
  %   sd=nanstd(d);
  %   finIx=isfinite(d);
  %   n=sum(finIx);
  %   sem=sd./sqrt(n);
  %   ci=sem.* -permute(tinv(.05/2,permute(n,[3 2 1])-1),[3 2 1]);
  %   d_var=ci;
  %   % d_var=sem;
  % end
  
  % concatenate, pick column (indep par) of interest, swapping columns and slices
  collCatD=cat(1,collCatD,permute(d(:,anColumnIx,:),[1 3 2]));
  
end


% **** jobs proper ****
switch job
  case 'PCA'
    [pc,transD,varnc]=PCexplore(collCatD,'tag',tag,'nPC',nPC,'plotType',{'scatter','density'});
    disp('parameters contributing to PC (most important for PC1 first):'),
    [nix,ix]=sort(abs(pc(:,1)),1,'descend');
    disp([strvcat(par(ix)) repmat('  ',numel(ix),1) num2str(pc(ix,1:3))]);
  otherwise 
    error('illegal job')
end


% ----- printing -----
if ~isempty(printas), 
  print(printas,figName,'-r400'); 
end