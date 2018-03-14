% This is the basic template for a script which collects, analyzes and
% plots processed FIELD POTENTIAL data from in vitro drug experiments. In
% PART I the type of analysis is set up. You have to define independent
% variables (e.g. concentration) and dependent variables (the things we
% measure, e.g. firing rates). Furthermore, and importantly, you may (have
% to) select specific parts of the data in the data base. This is what
% 'pa.subsetPar' is good for. Details see below. PART II is the specific
% analysis.

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
  };

% DEPENDENT PARAMETERS TO BE PLOTTED
% name of parameter, normalize?
% (abbreviations: std=standard deviation, cv=coefficient of variation)
ap.depPar={...
  % ------- FP-related parameters -------------
  'autoCDecay',           true;... % decay of autocorrelation of envelope of raw data
  'mnBurstPeak',          true;... % mean peak of envelope of burst (within first 500 ms)
  'mdBurstPeak',          true;... % median peak of envelope of burst (within first 500 ms)
  'cvBurstPeak',          false;... % cv of peak of envelope of burst (within first 500 ms)
  'asBurstPeak',          false;... % asymmetry of peak of envelope of burst (within first 500 ms)
  'mnBurstStrength',      true;... % mean burst integral (=area above base line; unit: mV*s)
  'mdBurstStrength',      true;... % median burst integral
  'cvBurstStrength',      false;... % cv of burst integral
  'asBurstStrength',      false;... % asymmetry of burst integral
  'burstIntegral',        true;... % sum of burst integrals divided by recording time (unit: mV)
  'relTimeInBurst',       true;... % relative time spent in bursts ('active' time)
  'mnBurstLen',           true;... % mean length of bursts
  'mdBurstLen',           true;... % median of length of bursts
  'maxBurstLen',          true;... % length of longest burst
  'cvBurstLen',           false;... % cv of length of bursts
  'asBurstLen',           false;... % asymmetry of length of bursts: (95th percentile-median)/(median-5th percentile)
  'burstRate',            true;... % burst rate
  'mnSilentPerLen',       true;... % mean length of silent period
  'mdSilentPerLen',       true;... % median of length of silent period
  'maxSilentPerLen',      true;... % length of longest silent period
  'cvSilentPerLen',       false;... % cv of length of silent period
  'asSilentPerLen',       false;... % asymmetry of length of bursts: (95th percentile-median)/(median-5th percentile)
  'specpower',            true;... % power in freq band to be specified below
  % ------- SPX-related parameters -------------
  'eventRate',            true;... % firing rates
  'mnPethEarlyFr',        true;... % mean of early firing rate
  'cvPethEarlyFr',        false;... % cv of early firing rate
  'mnPethIntermedFr',     true;... % mean of intermediate firing rate
  'cvPethIntermedFr',     false;... % cv of intermediate firing rate
  'mnPethLateFr',         true;... % mean of late firing rate
  'cvPethLateFr',         false;... % cv of late firing rate
  };

% DATA DISPLAY OPTIONS
% if true, data for each individual experiment will be plotted
ap.IndividExperimentPlot=true;
% output format for figures (e.g. '-djpeg95', '-dpsc2'; set to [] if
% figures shall not be saved)
ap.printFig='-djpeg95';
ap.printFig=[];

% ANALYSIS PARAMETERS
% --- LFP SECTION
% i) data preprocessing 
% the frequency to which LFP data will be down- or upsampled prior to all
% analyses
ap.sampleFreq=500;
% lowpass cutoff frequency (Hz)
ap.loCFreq=100;
% highpass cutoff frequency (Hz)
ap.hiCFreq=.2;
% ii) autocorrelation parameters
% lag (ms)
ap.autoCLag=10000;
% iii) spectral analysis
% range of frequencies to be represented in power spectra
ap.psdFreq=[4 80];
% frequency band from which to compute power (parameter 'specpower' above)
ap.pFreq=[40 80];
% length (ms) beyond which a burst will be regarded 'long' and subjected to
% spectral analysis - may be set to 'original'
ap.minLongBurstLength=2500;
% length of elementary segment from which to compute psd (ms)
ap.ppSeg=1024;
% extension of burst cutout window (ms) relative to the beginning and end
% of each burst. For example, [-100 100] adds 100 ms pre-burst data and 100
% ms post-burst data. [0 0] takes each burst as it had been detected. 
% Applies to spectral analysis only.
ap.winBuCutout=[0 0]

% --- SPX SECTION
% iv) peri-event time histograms
% ms, the time interval around a burst start (t=0) in which spikes will be
% collected. ** NOTE: the interval and bin width (below) must be chosen
% such that a bin with exactly zero as its left border results 
ap.pethIntv=[-50 1000];
% ms, bin width of the peth
ap.pethBinW=5;
% ms, intervals defining the 'early', 'intermediate' and 'late' firing rates
ap.pethEarlyFrIntv=[0 40];
ap.pethIntermedFrIntv=[40 250];
ap.pethLateFrIntv=[500 1000];


% PATHS TO DATA DIRECTORIES
% - rootp points to the *base* directory containing results from spike &
% burst detection
ap.rootPath='c:\hh\_tbs\bu\';

% resPath points to the directory in which figures and final results file
% will be collected
ap.resPath=ap.rootPath;
% name of the *.mat file (without extension) into which results will be
% saved
ap.resFn='stuffX'

% the name of the 'master table' (excel) file. 
% ** only the first worksheet will be read **
ap.masterFn=[ap.rootPath '\stuffX.xls'];

% -----------------------------------------------------------------------
% ----- PART II: call the data-reading & -processing functions ----------
% -----------------------------------------------------------------------
% retrieve information about the columns (=independent parameters) in the
% *.xls file specified above and perform some checks
indepPar=tslbatch_defdata(ap);

% read data
indepPar=tslbatch_readwsdata(ap.masterFn,indepPar);

% % write raw data into mat files
% tslbatch_abf2mat(indepPar,ap);

% process data
tslbatch_procdata_fp(indepPar,ap);

% plot SPX-related results
tslbatch_plotdata([ap.resPath  ap.resFn]);

% plot FP-related results
tslbatch_plotdata_fp([ap.resPath  ap.resFn]);
