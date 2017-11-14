% set paths
compName=lower(getenv('computername'));
switch compName
  case {'hh-i7'}
    ds.dataPath='d:/hh/projects/ctx_propoSevo/data_extra/';
    ds.figPath='d:/hh/projects/ctx_propoSevo/rawFig/';
  case {'hh64'}
    ds.dataPath='d:/hh/projects/ctx_propoSevo/data_extra/';
    ds.figPath='d:/hh/projects/ctx_propoSevo/rawFig/';
  otherwise
    error('machine not defined');
end

% full path & name of mat file (without extension)
ds.mFnFull=[ds.dataPath  '/PropoSevo_ctx_wtb6_SPX_sum_norm'];

% levels of independent variable (=concentrations in most cases) to deal
% with IN ASCENDING ORDER
ds.curIpVal=[0 .2 .5 .7];


% the independent var value to which data will be normalized (set to [] for
% no normalization)
ds.normIpVal=0;
% ds.normIpVal=[];

% average electrodes from cultures 
ds.doAvEl=true;

% outliers can be eliminated in either of two ways: 
% - if outlierType is set to 'abs' outlierThresh MUST contain TWO values in
% ASCENDING order, namely the lower and upper limit of acceptable values
% - if outlierType is set to 'std' outlierThresh MUST contain ONE value,
% the logic being that for each independent value (concentration) data
% points beyond [mean +/- outlierThresh*standard deviation of the whole population]
% will be kicked
% In case you don't want to eliminate any data points set outlierType to
% 'none'
ds.outlierType='std';
ds.outlierThresh=2.5;

% ds.outlierType='abs';
% ds.outlierThresh=[0 3];
% 
ds.outlierType='none';

% parameter(s) of interest
% ds.par={'mdSilentPerLen'};
% ds.par={'maxSilentPerLen'};
% ds.par={'mdBurstLen'};
% ds.par={'asBurstLen'};
% ds.par={'maxBurstLen'};

% ds.par={'burstRate'};
% ds.par={'burstIntegral'};

% ds.par={'specpower'};
% ds.par={'mnBurstStrength_norm'};
% ds.par={'mnBurstPeak'};
% ds.par={'relTimeInBurst'};

% ds.par={'eventRate_norm'};
ds.par={'eventRate'}; 
% ds.par={'cvISI'}; 
% ds.par={'eventRateInBurst'};
% ds.par={'fractionEvOutBurst'};
% ds.par={'fractionEvInBurst'};
% ds.par={'mnPethEarlyFr'};
% ds.par={'mnPethIntermedFr'};
% ds.par={'mnPethLateFr'};

% ds.par={'mdSTA'};


ds.printas=[];  % don't save graphics
% ds.printas='-djpeg97'; % save as jpeg ('working style')
% ds.printas='-dpsc2'; % save as color post script (nice finish)

% **** job to be performed & specific options ***
ds.job='list';
ds.job='dose response graph';
% ds.job='effect dependence';
% ds.job='fit';

ds.jobOpt={''}
switch ds.job
  
  case 'fit';
    ds.jobOpt={'Hill'};
    ds.jobOpt={'Hillbilly'};
    ds.jobOpt={'linear'};

  case 'dose response graph'
    % create an unfitted and/or 2D dose response graph
    % By default, unmatched samples will be plotted in a line graph. If a
    % '.._norm' parameter is chosen the data are by definition matched. For
    % options, use any combination of the options
    %       {'matched samples','bar','list'};
    % ** note: if two parameters are specified a 2D plot will be created
    % this will kick any recording in which there is no data for any of the
    % independent values
    ds.jobOpt={'matched samples'};
    ds.jobOpt={''};
    % ds.jobOpt={'list'};
    
  case 'effect dependence'
    % the parameter of interest (variable 'ds.par' defined above) will be
    % plotted versus some property of the recording. The 'property' may be
    % any of the precomputed dependent parameters or something else (see
    % below) at a (drug) condition of your choice.
    % Choose either a precomputed dependent parameter like 'mnBurstLen' or
    % 'burstiness'
    ds.jobOpt={'mdBurstLen'};
    % the (drug) condition of the property in question
    ds.curAbscissIpVal=0;

  case 'list'

  otherwise
    error('illegal job')
end

% graphics options
labelscale('fontSz',12,'scaleFac',.8,'lineW',2,'markSz',9);

% ** call datasniffer (requesting all output args)
[d,expName,indepParLevel]=datasniffer_tslbatch_01(ds);