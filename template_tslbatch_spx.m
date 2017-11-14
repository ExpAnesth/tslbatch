% This is the basic template for a script which collects, analyzes and
% plots processed SPIKE data from in vitro drug experiments. In PART I the
% type of analysis is set up. You have to define independent variables
% (e.g. concentration) and dependent variables (the things we measure, e.g.
% firing rates). Furthermore, and importantly, you may (have to) select
% specific parts of the data in the data base. This is what 'ap.subsetPar'
% is good for. Details see below. PART II is the specific analysis.

% -----------------------------------------------------------------------
% ----- PART I: setup of data information & analysis parameters ---------
% -----------------------------------------------------------------------

% INDEPENDENT PARAMETER OF INTEREST (ONE ONLY)
% name of parameter, control value (i.e. the value whose firing rates etc.
% will be used for normalization; usually this is zero)
ap.indepPar={...
  'drug1_conc', 0;...
  };

% INDEPENDENT PARAMETERS RESTRICTING THE SUBSET OF DATA TO BE ANALYZED
% name of parameter, min, max, specific values to be purged
ap.subsetPar={...
  'drug1_conc', 0, 0.2, [-1];...
  'drug1_applicRank',0, 2, [-1];...
  'undef1',0, 1, [-1];...
  };

% DEPENDENT PARAMETERS TO BE PLOTTED
% name of parameter, normalize?
% (abbreviations: std=standard deviation, cv=coefficient of variation)
ap.depPar={...
  'eventRate',            true;... % firing rates
  'cvISI',                false;...% cv of inter-spike-interval (makes sense only for single units)
  'fractionEvInBurst',    true;... % fraction of spx in bursts
  'relTimeInBurst',       true;... % relative time spent in bursts ('active' time)
  'mnNEvPerBurst',        true;... % mean of number of events in burst
  'cvNEvPerBurst',        true;... % cv of number of events in burst
  'mnIntraBurstEvRate',   true;... % mean of intra-burst event rate
  'cvIntraBurstEvRate',   true;... % cv of intra-burst event rate
  'mnBurstLen',           true;... % mean length of bursts
  'mdBurstLen',           true;... % median of length of bursts
  'cvBurstLen',           true;... % cv of length of bursts
  'asBurstLen',           true;... % asymmetry of length of bursts: (95th percentile-2*median+5th percentile)/median
  'burstRate',            true;... % burst rate
  'mnSilentPerLen',       true;... % mean length of silent period
  'mdSilentPerLen',       true;... % median of length of silent period
  'cvSilentPerLen',       true;... % cv of length of silent period
  'asSilentPerLen',       true;... % asymmetry of length of silent periods
  'mnPethEarlyFr',        true;... % mean of early firing rate
  'cvPethEarlyFr',        true;... % cv of early firing rate
  'mnPethIntermedFr',     true;... % mean of intermediate firing rate
  'cvPethIntermedFr',     true;... % cv of intermediate firing rate
  'mnPethLateFr',         true;... % mean of late firing rate  
  'cvPethLateFr',         true;... % cv of late firing rate    
};


% DATA DISPLAY OPTIONS
% if true, data for each individual experiment will be plotted
ap.IndividExperimentPlot=true;
% output format for figures (e.g. '-djpeg95', '-dpsc2'; set to [] if
% figures shall not be saved)
ap.printFig='-djpeg95';
% ap.printFig=[];

% ANALYSIS PARAMETERS
% i) general burst detection parameters
% ms, the minimal length of a silent period before a burst
% - may be set to 'original'
ap.minSilentPerDur='original';
% maximal number of events in pre-burst silent period
% - may be set to 'original'
ap.minSilentPerNEv='original';
% length (ms) beyond which a burst will be regarded a 'freak' event 
% - may be set to 'original'
% ** NOT YET IMPLEMENTED **
ap.maxNormalBurstLength=inf;


% ii) peri-event time histograms
% ms, the time interval around a burst start (t=0) in which spikes will be
% collected. ** NOTE: the interval and bin width (below) must be chosen
% such that a bin with exactly zero as its left border results 
ap.pethIntv=[-50 1100];
% ms, bin width of the peth
ap.pethBinW=5;
% ms, intervals defining the 'early', 'intermediate' and 'late' firing rates
ap.pethEarlyFrIntv=[0 50];
ap.pethIntermedFrIntv=[50 200];
ap.pethLateFrIntv=[500 1000];


% PATHS TO DATA DIRECTORIES
% - rootp points to the *base* directory containing results from spike &
% burst detection
ap.rootPath='d:\_data\otc_ctx\bicuIncub\';

% resPath points to the directory in which figures and final results file
% will be collected
ap.resPath=ap.rootPath;
% name of the *.mat file (without extension) into which results will be
% saved
ap.resFn='bicuIncub_Mg'

% the name of the 'master table' (excel) file. 
% ** it must reside in the  root directory as defined above **
% ** only the first worksheet will be read **
ap.masterFn=[ap.rootPath '\masterTableBicuMg.xls'];

% -----------------------------------------------------------------------
% ----- PART II: call the data-reading & -processing functions ----------
% -----------------------------------------------------------------------
% retrieve information about the columns (=independent parameters) in the
% *.xls file specified above and perform some checks
indepPar=tslbatch_defdata(ap);
% read data
indepPar=tslbatch_readwsdata(ap.masterFn,indepPar);
% process data 
tslbatch_procdata(indepPar,ap);
% plot them
tslbatch_plotdata([ap.resPath  ap.resFn]);
