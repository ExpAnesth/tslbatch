function tslbatch_procdata_fp(indepPar,ap)
% ** function tslbatch_procdata_fp(indepPar,ap)
% core processing function of tslbatch workhorse.

% to do:
% differentiate between different values of 'isvalid' and perform analyses
% accordingly

% -------------------------------------------------------------------------
%             PART 0: preliminaries: auxil pars & checks of ap
% -------------------------------------------------------------------------

wp.unitList={...
  'µV', .001;
  'mV', 1;
  'V', 1000;
  'pA', 1;
  'nA',1000};
wp.scaleFac=1;
% internal code for long bursts
wp.isLongBurst=3;
% internal code determining whether FP, SPX or both shall be used
wp.useFPSPX=1;
wp.useFP=2;
wp.useSPX=3;
% histogram bins for envelope of field potential (needed for subtraction of
% base line in autocorr computation): noise is rarely above 50 uV peak to
% peak, and envelope is strictly positive
wp.envAmpBin=0:.0002:.03;
% interval in which to find peak of envelope of burst
wp.buPeakFindIntv=[-5 500];


etslconst;

% names & number of parameters in indepPar
parNm={indepPar.name};
nPar=length(indepPar);

% some indices into indepPar needed throughout
expIDIx=strmatch('expID',parNm);
isValidIx=strmatch('isValid',parNm);
fileNameIx=strmatch('fileName',parNm);
chanNameIx=strmatch('chanName',parNm);
drug1_concIx=strmatch('drug1_conc',parNm);
drug1_applicRankIx=strmatch('drug1_applicRank',parNm);
undef1Ix=strmatch('undef1',parNm);
undef2Ix=strmatch('undef2',parNm);

etslconst;

% - check ap & set flags
doUseOrig_minLongBurstLength=false;
doUseMeanSlopeThresh=false;
% this is a tribute to the fact that this parameter has been introduced
% quite late
if ~isfield(ap,'minLongBurstLength')
  ap.minLongBurstLength=inf;
end
% catch obvious blunders in intervals for spec analysis
if ap.minLongBurstLength<ap.ppSeg+diff(ap.winBuCutout)
  error(['minimal ''long'' burst length is too small']);
end

% - simple checks of peth parameters:
% check whether there's a bin with zero as its left border (while doing so
% generate a fake results array as a template)
[hstTemplate,pethBin]=tslpeth(1,1,'interval',ap.pethIntv,'binw',ap.pethBinW);
% ** hstTemplate, and thus the preallocated variable holding all
% histograms, must be all nans (as opposed to zeros) because this is how
% recordings with zero events should be represented
hstTemplate=nan*hstTemplate;
if isempty(find(pethBin==0))
  error('the bins for the PETH must include a bin with zero as its left border');
end
% see whether pefr intervals are included in pethBin and, while doing so,
% generate indices into pethBin
tmp={'ap.pethEarlyFrIntv','ap.pethIntermedFrIntv','ap.pethLateFrIntv'};
for g=1:3
  intv=eval(tmp{g});
  if intv(1)<ap.pethIntv(1) || intv(2)>ap.pethIntv(2)
    error([tmp{g}  ' must be completely enclosed by ap.pethIntv']);
  end
  pethBinIx{g}=find(pethBin>=intv(1) & pethBin<intv(2));
end

% make sure that all dependent parameters requested are part of the
% staple computed here. 
% The second column of masterDepPar indicates whether the dependent
% variable is scalar (1), or sth else (currently not used); the third
% column indicates whether the parameter is a FP-related (1) or a
% SPX-related one (2) - this information will be used later for the
% generation of summary plots

masterDepPar={...
  'autoCDecay',           1, 1;... % decay of autocorrelation
  'mnBurstPeak',          1, 1;... % mean peak of envelope of burst (within wp.buPeakFindIntv ms)
  'mdBurstPeak',          1, 1;... % median peak of envelope of burst (within wp.buPeakFindIntv ms)
  'cvBurstPeak',          1, 1;... % cv of peak of envelope of burst (within wp.buPeakFindIntv ms)
  'asBurstPeak',          1, 1;... % asymmetry of peak of envelope of burst (within wp.buPeakFindIntv ms)
  'mnBurstStrength',      1, 1;... % mean burst integral (=area above base line; unit: mV*s)
  'mdBurstStrength',      1, 1;... % median burst integral
  'cvBurstStrength',      1, 1;... % cv of burst integral
  'asBurstStrength',      1, 1;... % asymmetry of burst integral  
  'burstIntegral',        1, 1;... % sum of burst integrals divided by recording time (unit: mV)  
  'relTimeInBurst',       1, 1;... % relative time spent in bursts ('active' time)
  'mnBurstLen',           1, 1;... % mean length of bursts
  'mdBurstLen',           1, 1;... % median of length of bursts
  'maxBurstLen',          1, 1;... % length of longest burst
  'cvBurstLen',           1, 1;... % cv of length of bursts
  'asBurstLen',           1, 1;... % asymmetry of length of bursts: (95th percentile-median)/(median-5th percentile)
  'burstRate',            1, 1;... % burst rate
  'mnSilentPerLen',       1, 1;... % mean length of silent period
  'mdSilentPerLen',       1, 1;... % median of length of silent period
  'maxSilentPerLen',      1, 1;... % length of longest silent period  
  'cvSilentPerLen',       1, 1;... % cv of length of silent period
  'asSilentPerLen',       1, 1;... % asymmetry of length of silent period
  'specpower',            1, 1;... % power in freq band to be specified below
  'eventRate'           , 1, 2;... % firing rates
  'mnPethEarlyFr'       , 1, 2;... % mean of early firing rate
  'stdPethEarlyFr'      , 1, 2;... % std of early firing rate
  'cvPethEarlyFr'       , 1, 2;... % cv of early firing rate
  'mnPethIntermedFr'    , 1, 2;... % mean of intermediate firing rate
  'stdPethIntermedFr'   , 1, 2;... % std of intermediate firing rate  
  'cvPethIntermedFr'    , 1, 2;... % cv of intermediate firing rate
  'mnPethLateFr'        , 1, 2;... % mean of late firing rate  
  'stdPethLateFr'       , 1, 2;... % std of late firing rate
  'cvPethLateFr'        , 1, 2;... % cv of late firing rate    
};

tmp=setdiff(ap.depPar(:,1),masterDepPar(:,1));
if ~isempty(tmp)
  errordlg({'Analysis parameters'; ' '; strvcat(tmp); ' '; 'do not exist (check spelling - NO BLANKS ARE ALLOWED).'; ' '; 'Here is the full list of legal parameters:'; ' '; strvcat(masterDepPar(:,1))});
  error('see error dialog');
end
nMasterDepPar=size(masterDepPar,1);
nDepPar=size(ap.depPar,1);

% do autocorr, spectral analysis, burst analysis requiring raw data and spx
% analysis depending on choice of plot parameters (the standard burst
% analysis will be done unconditionally) 
doAutoCorrAnalysis=ismember('autoCDecay',ap.depPar(:,1));
doSpecAnalysis=ismember('specpower',ap.depPar(:,1));
doBurstIntegralAnalysis=any(ismember(...
  {'mnBurstPeak',...
  'mdBurstPeak',... 
  'cvBurstPeak',... 
  'asBurstPeak',... 
  'mnBurstStrength',...
  'mdBurstStrength',... 
  'cvBurstStrength',... 
  'asBurstStrength',... 
  'burstIntegral'},... 
  ap.depPar(:,1)));
doSpxAnalysis=any(ismember(...
  {'eventRate'        ,... 
  'mnPethEarlyFr'     ,... 
  'stdPethEarlyFr'    ,... 
  'cvPethEarlyFr'     ,... 
  'mnPethIntermedFr'  ,... 
  'stdPethIntermedFr' ,... 
  'cvPethIntermedFr'  ,... 
  'mnPethLateFr'      ,... 
  'stdPethLateFr'     ,... 
  'cvPethLateFr'}     ,... 
  ap.depPar(:,1)));

% -------------------------------------------------------------------------
%                PART I: purges of observations
% -------------------------------------------------------------------------
% ** the information about exclusion criteria has been transferred from
% ap.subsetPar to indepPar at this point (in tslbatch_defdata) **
% logical array pointing to observations to be killed
purgeArr=zeros(size(indepPar(1).d));
for g=1:length(indepPar)
  aPar=indepPar(g).d;
  if isnumeric(aPar)
    if strcmp(indepPar(g).dType,'char')
      error(['independent paramater ' indepPar(g).name ' is numeric, but should be char']);
    end
    aPar=double(indepPar(g).d);
    if isfinite(indepPar(g).min)
      purgeArr=purgeArr | aPar<indepPar(g).min;
    end
    if isfinite(indepPar(g).max)
      purgeArr=purgeArr | aPar>indepPar(g).max;
    end
    % if .purgeVal contains floats we have to take some precaution
    % against glitches due to roundoff errors/limited number precision.
    % Provisional solution here: take up to fourth decimal digit into account 
    if ~isempty(indepPar(g).purgeVal)
      purgeArr=purgeArr | ismember(round(aPar*1000),round(indepPar(g).purgeVal*1000));
    end
  else
    % if parameter is non-numeric 
    % - .max and .min are meaningless so no need to check them
    % - we can use ismember right away to identify undesired entries
    if ~isempty(indepPar(g).purgeVal)
      purgeArr=purgeArr | ismember(aPar,indepPar(g).purgeVal);
    end
  end
end

tmp=numel(find(~indepPar(isValidIx).d));
disp(['** excluding ' int2str(tmp) ' observations from analysis according to entries in worksheet master list']);
purgeArr=purgeArr | ~indepPar(isValidIx).d; 
disp(['** excluding a further ' int2str(numel(find(purgeArr))-tmp) ' observations from analysis according to entries in ap.subsetPar'])
% ** the great purge
if any(purgeArr)
  for g=1:length(indepPar)
    indepPar(g).d(purgeArr)=[];
  end
end

% -------------------------------------------------------------------------
%                PART II: setup of master results variables
% -------------------------------------------------------------------------
% identify individual 'experiments' 
[uExp,ix1,ix2]=unique(cellstr(indepPar(expIDIx).d));
nExp=length(uExp);

% some additional precautionary checks
for g=1:nExp
  if numel(unique(indepPar(chanNameIx).d(ix2==g)))>2
    warndlg(['there are more than 2 channels listed for experiment ' uExp(g) ' - you may want to make sure this is correct']);
  end
end

% identify the independent parameters to be used for analysis
aIndepIx=find([indepPar.use]);
nIndepPar=length(aIndepIx);
if nIndepPar>1
  error([mfilename ' works only with one independent parameter (check ''ap.indepPar'')']);
end

% *** if this is an experiment with drug1_conc as the indep var...
if strcmp(indepPar(aIndepIx).name,'drug1_conc')
  % set drug 1 applic rank order of wash to 0
  tmpIx=indepPar(drug1_concIx).d==-1;
  indepPar(drug1_applicRankIx).d(tmpIx)=0;
  % check whether drug 1 applic rank order of control (as defined by conc==0) is 0
  tmpIx=indepPar(drug1_concIx).d==0 & indepPar(drug1_applicRankIx).d~=0;
  if any(tmpIx)
    warndlg({'In ';' '; strvcat(indepPar(expIDIx).d{tmpIx}); ' '; 'at least one control recording (defined as [drug1]=0) has an application rank different from 0 - make sure that this is correct';'(If so and you don''t want this message to keep popping up talk to HH'});
  end
end

% which and how many levels present in data?
indepParLevel=unique(indepPar(aIndepIx).d);
indepParNLevel=numel(indepParLevel);
% which of these levels is the one to be used for normalization?
indepParNormIx=find(indepParLevel==indepPar(aIndepIx).normVal);

% number of observations (=atsl files)
nObs=size(indepPar(1).d,1);
% convert channel name to ID...
[nix,nix2,tmpChanID]=unique(indepPar(chanNameIx).d);
% ...to obtain number of unique combinations of experiment ID and channel
nUExpChan=size(unique([ix2 tmpChanID],'rows'),1);

% put out some information
disp(['** independent parameter: ' indepPar(aIndepIx).name '; levels: ' num2str(indepParLevel(:)')])
disp(['** ' int2str(nExp) ' experiments, ' int2str(nObs) ' observations'])


% ---------- specific preparatory works for analyses --------------
% mock spectral results pars
[psdPTemplate,psdF]=fspecp(rand(cont2discrete(ap.ppSeg,1000/ap.sampleFreq)+1,1),...
  1e6/ap.sampleFreq,'meth','fft','win',ap.ppSeg,'limFreq',ap.psdFreq);
psdPTemplate=psdPTemplate*nan;
% index to frequency band of interest
freqIx=psdF>=ap.pFreq(1) & psdF<=ap.pFreq(2);
% scaling factor for computation of power in specific band from psd
powScaleFac=ap.sampleFreq/(2*round(ap.ppSeg/ap.sampleFreq*1000));

% ---------- preallocation of processed results structure --------------
% --- dimensions of cell array containing variably sized values (tsl):
%   <one row> | indep par value | <unique combin experiment, channel>
vsTemplate=cell([1 indepParNLevel nUExpChan]);
% it is fortunate for various reasons if the individual cells of
% tsl-containing variables contain zero-row templates of etsl
[vsTemplate{:}]=deal(zeros(0,etslc.nCol));
% note UPPERCASE:
% - FP bursts
Etsl=vsTemplate;
% - FP silent periods
SilentEtsl=vsTemplate;
% - SPX 
Tsl=vsTemplate;

% --- dimensions of autocorrelation:
autoCLag_pts=cont2discrete(ap.autoCLag,1000/ap.sampleFreq);
%  lag (one-sided) | indep par value | <unique combin experiment, channel>
autoC=repmat(nan,[autoCLag_pts+1 indepParNLevel nUExpChan]);

% --- dimensions of psd variables:
%  freq bins | indep par value | <unique combin experiment, channel>
psDensMn=repmat(psdPTemplate,[1 indepParNLevel nUExpChan]);
psDensStd=psDensMn;

% --- dimensions of histogram variables:
%   peri-event time bins| indep par value | <unique combin experiment, channel>
pethTemplate=repmat(hstTemplate,[1 indepParNLevel nUExpChan]);
pethMn=pethTemplate;
pethStd=pethTemplate;

% --- dimensions of SCALAR values like burst length:
%   <one row> | indep par value | <unique combin experiment, channel>
scalarTemplate=repmat(nan,[1 indepParNLevel nUExpChan]);

% - decay of autocorr
autoCDecay=scalarTemplate;
% - relative time spent in bursts ('active' time)
relTimeInBurst=scalarTemplate;

% - mean, std, cv, etc. of peak of envelope of burst (within wp.buPeakFindIntv ms; unit: mV*s)
mnBurstPeak=scalarTemplate;
mdBurstPeak=scalarTemplate;     
stdBurstPeak=scalarTemplate;     
cvBurstPeak=scalarTemplate;
asBurstPeak=scalarTemplate;

% - mean, std, cv, etc. of burst integral (=area above base line; unit: mV*s)
mnBurstStrength=scalarTemplate;
mdBurstStrength=scalarTemplate;  
stdBurstStrength=scalarTemplate;     
cvBurstStrength=scalarTemplate;
asBurstStrength=scalarTemplate;
% sum of burst integrals divided by recording time (unit: mV)  
burstIntegral=scalarTemplate;        

% - mean, std, cv, etc. of length of bursts
mnBurstLen=scalarTemplate;
mdBurstLen=scalarTemplate;
maxBurstLen=scalarTemplate;
stdBurstLen=scalarTemplate;
cvBurstLen=scalarTemplate;
asBurstLen=scalarTemplate;
% - burst rate
burstRate=scalarTemplate;
% - mean, std, cv, etc. of length of silent period
mnSilentPerLen=scalarTemplate;
mdSilentPerLen=scalarTemplate;
maxSilentPerLen=scalarTemplate;
stdSilentPerLen=scalarTemplate;
cvSilentPerLen=scalarTemplate;
asSilentPerLen=scalarTemplate;
specpower=scalarTemplate;

% - firing rates
eventRate=scalarTemplate;
% - mean, std and cv of early, intermediate and late fr
mnPethEarlyFr=scalarTemplate;
stdPethEarlyFr=scalarTemplate;
cvPethEarlyFr=scalarTemplate;
mnPethIntermedFr=scalarTemplate;
stdPethIntermedFr=scalarTemplate;
cvPethIntermedFr=scalarTemplate;
mnPethLateFr=scalarTemplate;
stdPethLateFr=scalarTemplate;
cvPethLateFr=scalarTemplate;

% also export application rank order, undef1 and undef2 (may be useful)
drug1_applicRank=scalarTemplate;
undef1=scalarTemplate;
undef2=scalarTemplate;
% important: have a 1D cell array with as many elements as the results
% variables above have 'slices': it contains strings identifying the
% nUExpChan 'experiments'
expChanName=cell(1,nUExpChan);
% we also need rudimentary file info
fileInfo.dDir='';
fileInfo.fileName='';
fileInfo.chName='';
fileInfo.recTime=nan;
fileInfo.si=nan;
fileInfo.winBuCutout=ap.winBuCutout;
fileInfo=repmat(fileInfo,[1 indepParNLevel nUExpChan]);

% internal control: check whether all dependent variables are present at
% this point
s=whos;
tmp=setdiff(masterDepPar(:,1),{s.name});
if ~isempty(tmp)
  error('internal: preallocation failure - tell the programmer');
end

% -------------------------------------------------------------------------
%                PART III: collect data
% -------------------------------------------------------------------------
% the 'slice' index into the master* variables
masterSliceIx=0;
% container for files not found:
% - main results file containing etsl for FP
filesNotFound.fpEtsl=cell(0);
% - either of the raw data files (mat or abf)
filesNotFound.raw=cell(0);
% - file containing *tsl* from spx
filesNotFound.spxTsl=cell(0);
% - §§§ file containing *etsl* from spx
filesNotFound.spxEtsl=cell(0);
% container for files with empty FP/MP etsl 
filesEmptyEtsl=cell(0);
% container for files with empty SPX tsl 
filesEmptyTsl=cell(0);
for g=1:nExp
  % index into indepPar(expIDIx).d for current experiment
  curExpIx=find(ix2==g);
  % experimental ID
  tmpStr=indepPar(expIDIx).d{curExpIx(1)};
  % subdirectory in which data reside: 
  % i. date part of expID, underscore
  % ii. first three characters of expID identifying the setup
  datDir=[ap.rootPath tmpStr([6:16 1:4])];
  % unique channels of current experiment
  uChan=unique(indepPar(chanNameIx).d(curExpIx));
  nChan=numel(uChan); 
  expChan=[];
  % --- loop over channels
  for h=1:nChan
    % index into indepPar.d for channels of current experiment
    curExpChanIx=intersect(curExpIx,strmatch(uChan(h),indepPar(chanNameIx).d,'exact'));
    % current deblanked channel name
    curDeblankChName=uChan{h};
    curDeblankChName=curDeblankChName(~isspace(curDeblankChName));
    masterSliceIx=masterSliceIx+1;
    % fill expChanName: experiment ID + file name
    expChanName{masterSliceIx}=[indepPar(expIDIx).d{curExpChanIx(1)},...
      ', ', indepPar(chanNameIx).d{curExpChanIx(1)}];
    % --- for given channel, loop over files (=levels of indep par)
    for i=1:numel(curExpChanIx)
      % §§ number coding for usage of SPX, FP or both
      isValid=indepPar(isValidIx).d(curExpChanIx(i));
      % value of level
      curLevel=indepPar(aIndepIx).d(curExpChanIx(i));
      % index to column
      colIx=curLevel==indepParLevel;
      % set flags for the different analyses
      doAC=doAutoCorrAnalysis;
      doSA=doSpecAnalysis;
      doBIA=doBurstIntegralAnalysis;
      doSPX=doSpxAnalysis;
      % file name without extension as listed in master table
      curFnBase=indepPar(fileNameIx).d{curExpChanIx(i)};
      % name of FP results file
      curFnFpRes=[curFnBase '_' curDeblankChName '_FP_res'];
      % name of MP results file
      curFnMpRes=[curFnBase '_' curDeblankChName '_MP_res'];
      % name of SPX results file
      curFnSpxRes=[curFnBase '_' curDeblankChName '_SPX_res'];
      % name of cutouts file  (currently not needed)    
      % curFnFpCutout=[curFnBase '_' curDeblankChName '_FP_buCutout'];
      disp(['*** ' curFnBase '_' curDeblankChName '...']);

      % collect file name without extension & other info no matter whether
      % the data were found or not
      fileInfo(1,colIx,masterSliceIx).fileName=curFnBase;
      fileInfo(1,colIx,masterSliceIx).dDir=datDir;
      fileInfo(1,colIx,masterSliceIx).chName=curDeblankChName;      
      % also, collect all independent parameters to be saved in a format
      % identical to the scalar results variables
      drug1_applicRank(1,colIx,masterSliceIx)=indepPar(drug1_applicRankIx).d(curExpChanIx(i));
      undef1(1,colIx,masterSliceIx)=indepPar(undef1Ix).d(curExpChanIx(i));
      undef2(1,colIx,masterSliceIx)=indepPar(undef2Ix).d(curExpChanIx(i));      
      % -------------------------------------------------------------------
      % ---------------------- reading data -------------------------------
      % -------------------------------------------------------------------
      
      % ------- load FP/MP burst results & general information ------------
      fExist=[exist([datDir '\' curFnFpRes '.mat'],'file'),...
        exist([datDir '\' curFnMpRes '.mat'],'file')];
      % *** by default assume that the data are from extracellular
      % recordings (§ parameter not yet used)
      isExtraRec=true;
      if any(fExist)
        if fExist(1)
          load([datDir '\' curFnFpRes '.mat']);
        else
          load([datDir '\' curFnMpRes '.mat']);
          isExtraRec=false;
          doBIA=false;
          doAC=false;
          doSA=false;
        end          
        % ** collect vital information on recording: **
        % recording time in s
        fileInfo(1,colIx,masterSliceIx).recTime=head.ds.fileInfo.recTime;
        % §§§ sampling interval of downsampled (interpolated) data used here
        fileInfo(1,colIx,masterSliceIx).si=1e6/ap.sampleFreq;
        % detach etsl, silentEtsl and stats
        etsl=bu.etsl;
        silentEtsl=bu.silentEtsl;
        stats=bu.stats;
        if isempty(etsl)
          warning('no FP or MP bursts detected');
          filesEmptyEtsl{end+1}=curFnBase;
        end
        if ~isfield(stats,'cIntv')
          % if this is not an empty etsl...
          if ~isempty(etsl)
            % infer cIntvSpan from relative time spent in active state
            % (imprecise)
            cIntvSpan=sum(etsl(:,etslc.durCol))/stats.relTimeInBurst;
          else
            % empty etsl means no burst detected (not necessarily an
            % invalid recording), so don't put an error or warning
            cIntvSpan=nan; % diff(head.ds.fileInfo.recTime)*1000;
          end
        else
          cIntvSpan=diff(stats.cIntv);
        end
      else
        % §§ there's an ambiguity here: it cannot be resolved whether the
        % file looked for but not found is a FP or a MP one. Collect them
        % in the list of missing FP files for now.
        filesNotFound.fpEtsl{end+1}=curFnFpRes;
        etsl=[];
        silentEtsl=[];
      end
      
      % at this point heed the value of isValid
      switch isValid
        case 1
          % fp and spx analysis
        case 2
          % only fp analysis
          doSPX=false;
        case 3
          % only spx analysis
          doBIA=false;
          doAC=false;
          doSA=false;
          etsl=[];
          silentEtsl=false;
        otherwise
          error('internal: invalid value for ''use this recording''')
      end
      
      % -------------- raw or preconditioned FP data ----------------------
      if doBIA || doAC || doSA
        isExistentFile=[exist([datDir '\' curFnBase '.mat'],'file'),...
          exist([datDir '\' curFnBase '.abf'],'file')];
        if ~any(isExistentFile)
          filesNotFound.raw{end+1}=curFnBase;
          doBIA=false;
          doAC=false;
          doSA=false;
          % no need to do anything to collected results variables because
          % they are preallocated with nans (all other)
        else
          % ** file loading logics:
          % - if mat raw file exists load it
          % - if only abf file exists load it and preprocess data (time-
          % consuming)
          if isExistentFile(1)
            [d,si,fi]=matDload([datDir '\' curFnBase '.mat'],'channels',{curDeblankChName});
            % check whether filter freq and sampling freq match
            if ~isalmost(fi.si,1e6/ap.sampleFreq,.01)
              error('sampling frequency of *.mat file containing raw data does not correspond to currently chosen ap.sampleFreq - either change ap.sampleFreq or recreate raw data files with current ap.sampleFreq');
            end
          else
            [d,osi,fi]=abfload([datDir '\' curFnBase '.abf'],'channels',uChan(h));
            % filter
            if isfinite(ap.loCFreq)
              d=lofi(d,osi,ap.loCFreq);
            end
            % resample to target freq
            nPt=size(d,1);
            d=interp1((1:nPt)',d,linspace(1,nPt,nPt*osi*(ap.sampleFreq/1e6))');
            fi.si=1e6/ap.sampleFreq;
            % hipass filter
            if isfinite(ap.hiCFreq)
              d=hifi(d,fi.si,ap.hiCFreq);
            end
          end
          % possibly scale data to mV
          chIx=strmatch(uChan(h),fi.recChNames,'exact');
          ix=strmatch(fi.recChUnits{chIx},wp.unitList(:,1),'exact');
          % if unit not found don't scale data and warn
          if isempty(ix)
            warndlg(['physical units on chosen channel unknown (' fi.recChUnits{chIx} ') - data will not be scaled']);
          else
            wp.scaleFac=wp.unitList{ix,2};
          end
          % scale
          d=d*wp.scaleFac;
          % *** envelope of d ***
          dEnv=abs(hilbert(d));
        end
      end

      % ------- load SPX results  ---------------------------------
      tsl=[];
      if doSPX
        if exist([datDir '\' curFnSpxRes '.mat'],'file')
          % load only evt
          load([datDir '\' curFnSpxRes '.mat'],'evt');
          % accomodate old versions of results data in which neither tsl
          % nor etsl were cells; also, make sure data are not from multiple
          % single (separated) units
          if iscell(evt.tsl)
            if numel(evt.tsl)>1
              error('cannot deal with multiple single units')
            end
            tsl=evt.tsl{1};
          else
            tsl=evt.tsl;
          end
          clear evt
          if isempty(tsl)
            filesEmptyTsl{end+1}=curFnBase;
          end
        else
          filesNotFound.spxTsl{end+1}=curFnSpxRes;
        end
      end

      % -------------------------------------------------------------------
      % ---------------------- computations -------------------------------
      % -------------------------------------------------------------------

      % ---------------- burst analysis ---------------------
      if isempty(etsl)
        % if etsl is empty skip all burst analysis. All results variables
        % are preallocated with nans, so there will be no funny gaps. We
        % don't want to insert an empty matrix where there is already an
        % array of NaNs, so leave Etsl alone.
        % By extension, spectral analysis won't happen, either, because it
        % is done on bursts only
        doSA=false;
      else
        % --- analyses of burst peak/strength/integral
        if doBIA
          % - number of bursts
          nBu=size(etsl,1);
          % - max length of excerpts in which to find peaks in pts
          maxNExcPts=diff(cont2discrete(wp.buPeakFindIntv,1000/ap.sampleFreq,'intv',1))+1;
          % - lengths of bursts in ticks
          buL=cont2discrete(etsl(:,etslc.durCol),1000/ap.sampleFreq,'intv',1);
          % - beginnings of bursts in ticks
          buB=cont2discrete(etsl(:,etslc.tsCol),1000/ap.sampleFreq,'intv',1);
          % - cutout times according to wp.buPeakFindIntv in ticks
          buC=cont2discrete(repmat(etsl(:,etslc.tsCol),1,2)+...
            repmat(wp.buPeakFindIntv,nBu,1),1000/ap.sampleFreq,'intv',1);
          % curb 
          buC(buC<1)=1;
          buC(:,2)=min(buC(:,2),buB+buL);
          excPts=diff(buC,1,2)+1;
          % preallocate
          tmpBuD=repmat(nan,maxNExcPts,nBu);
          % - cut out burst beginnings and place in tmpBuD
          for bIx=1:nBu
            tmpBuD(1:excPts(bIx),bIx)=dEnv(buC(bIx,1):buC(bIx,2));
          end
          % determine max peak
          r=evdeal(tmpBuD,'idx',{'minmaxpeak'});
          mnBurstPeak(1,colIx,masterSliceIx)=nanmean(r.maxPeak);
          stdBurstPeak(1,colIx,masterSliceIx)=nanstd(r.maxPeak);
          tmpx=prctile(r.maxPeak,[5 50 95]);
          mdBurstPeak(1,colIx,masterSliceIx)=tmpx(2);
          asBurstPeak(1,colIx,masterSliceIx)=(tmpx(3)-tmpx(2))/(tmpx(2)-tmpx(1));
          % - place into corresponding column in etsl
          etsl(:,etslc.amplCol)=r.maxPeak';
          % compute integral of envelope, heeding units (mV*s):
          % - cumulative sum of whole data trace
          dd=cumsum(dEnv);
          % - e voila parameter 'burst strength'
          buI=(dd(buB+buL)-dd(buB))/ap.sampleFreq;
          % - mean, std, median and 'asymmetry' of this 
          mnBurstStrength(1,colIx,masterSliceIx)=mean(buI);
          stdBurstStrength(1,colIx,masterSliceIx)=std(buI);
          tmpx=prctile(buI,[5 50 95]);
          mdBurstStrength(1,colIx,masterSliceIx)=tmpx(2);
          asBurstStrength(1,colIx,masterSliceIx)=(tmpx(3)-tmpx(2))/(tmpx(2)-tmpx(1));
          % - sum of burst integrals divided by recording time (unit: mV)
          burstIntegral(1,colIx,masterSliceIx)=sum(buI)/(cIntvSpan/1000);
        end
        % now copy (silent) etsls into respective containers
        Etsl{1,colIx,masterSliceIx}=etsl;
        SilentEtsl{1,colIx,masterSliceIx}=silentEtsl;
        % ----- burst statistics: 
        % - relative time spent in bursts ('active' time)
        relTimeInBurst(1,colIx,masterSliceIx)=stats.relTimeInBurst;
        % -  mean, std length of bursts
        mnBurstLen(1,colIx,masterSliceIx)=stats.mnBurstLen;
        stdBurstLen(1,colIx,masterSliceIx)=stats.stdBurstLen;
        maxBurstLen(1,colIx,masterSliceIx)=max(etsl(:,etslc.durCol));
        
        % - burst rate
        burstRate(1,colIx,masterSliceIx)=stats.burstRate;
        % - mean, std length of silent period
        mnSilentPerLen(1,colIx,masterSliceIx)=stats.mnSilentPerLen;
        stdSilentPerLen(1,colIx,masterSliceIx)=stats.stdSilentPerLen;
        maxSilentPerLen(1,colIx,masterSliceIx)=max(silentEtsl(:,etslc.durCol));
        % - add statistics not present in struct stats (CVs will be 
        % computed outside loop):
        tmpx=prctile(etsl(:,etslc.durCol),[5 50 95]);
        mdBurstLen(1,colIx,masterSliceIx)=tmpx(2);
        asBurstLen(1,colIx,masterSliceIx)=(tmpx(3)-tmpx(2))/(tmpx(2)-tmpx(1));
        if ~isempty(silentEtsl)
          tmpx=prctile(silentEtsl(:,etslc.durCol),[5 50 95]);
          mdSilentPerLen(1,colIx,masterSliceIx)=tmpx(2);
          asSilentPerLen(1,colIx,masterSliceIx)=(tmpx(3)-tmpx(2))/(tmpx(2)-tmpx(1));
        end
      end

      % ---------------- autocorrelation analysis -------------------------
      if doAC
%         % compute base line of envelope and subtract, then compute autocorr
%         % (instead of autocov)
%         tmpN=histc(dEnv,wp.envAmpBin);
%         % peak of hist
%         [nix,ix]=max(tmpN);
        % autocorr
        tmpC=xcov(dEnv,autoCLag_pts,'coeff');
        % keep one-sided
        tmpC=tmpC(autoCLag_pts+1:end);
        autoC(:,curLevel==indepParLevel,masterSliceIx)=tmpC;
        % compute autoCDecay:
        lim=.25;
        tmp=find(tmpC(1:end-1)>lim & tmpC(2:end)<lim,1);
        if ~isempty(tmp)
          autoCDecay(1,colIx,masterSliceIx)=tmp/(ap.sampleFreq/1000);
        end
      end

      % ---------------- burst-based spectral analysis --------------------
      if doSA
        % file-specific checks of cutout lengths etc.
        etsl(:,etslc.tagCol)=0;
        % all bursts long enough (taking ap.winBuCutout into consideration) 
        % & not too close to borders
        tmpIx=etsl(:,etslc.durCol)+diff(ap.winBuCutout)>=ap.minLongBurstLength;
        tmpIx=tmpIx & ...
          etsl(:,etslc.tsCol)+ap.winBuCutout(1)>0 & ...
          sum(etsl(:,[etslc.tsCol etslc.durCol]),2)+ap.winBuCutout(2)< diff(fileInfo(1,colIx,masterSliceIx).recTime)*1000;
        etsl(tmpIx,etslc.tagCol)=wp.isLongBurst;
        % index to long bursts
        buLongIx=find(etsl(:,etslc.tagCol)==wp.isLongBurst);
        nBuLong=numel(buLongIx);
      end
      if doSA && nBuLong
        bup=repmat(psdPTemplate,[1 nBuLong]);
        % - beginning of interval containing burst in ticks
        buB=cont2discrete(etsl(:,etslc.tsCol)+ap.winBuCutout(1),1000/ap.sampleFreq,'intv',1);
        % - same for end
        buE=cont2discrete(sum(etsl(:,[etslc.tsCol etslc.durCol]),2)+ap.winBuCutout(2),1000/ap.sampleFreq,'intv',1);
        for bix=1:nBuLong
          bud=d(buB(buLongIx(bix)):buE(buLongIx(bix)));
          nPt=numel(bud);
          [bup(:,bix),f]=fspecp(bud,1e6/ap.sampleFreq,'meth','fft','win',ap.ppSeg,'limFreq',ap.psdFreq);
        end
        bupMn=mean(bup,2);
        bupStd=std(bup,0,2);
        psDensMn(:,curLevel==indepParLevel,masterSliceIx)=bupMn;
        psDensStd(:,curLevel==indepParLevel,masterSliceIx)=bupStd;
        % power in specific freq band
        specpower(1,colIx,masterSliceIx)=sum(bupMn(freqIx))*powScaleFac;
      end

      % ---------------- spikes analysis ----------------------------------
      if doSPX && ~isempty(tsl)
        Tsl{1,colIx,masterSliceIx}=tsl;
        % *** in contrast to burst parameters compute average firing rate
        % as spike count divided by total recording time (NOT cIntvSpan)
        eventRate(1,colIx,masterSliceIx)=numel(tsl)/diff(head.ds.fileInfo.recTime);
        if ~isempty(etsl)
          % - peth: collect spx/bin, not spx frequency; do not collect spx
          % belonging to next ref event
          [hst,pethBin]=tslpeth(etsl(:,etslc.tsCol),tsl,'interval',ap.pethIntv,...
            'intvDistance',.1,'binw',ap.pethBinW,'convHz',0);
          % - within each experiment, average (over bursts) histograms
          % ** note that tslpeth spits out nans for bursts too close to
          % either border **
          hstMn=nanmean(hst,2);
          hstStd=nanstd(hst,0,2);
          pethMn(:,curLevel==indepParLevel,masterSliceIx)=hstMn;
          pethStd(:,curLevel==indepParLevel,masterSliceIx)=hstStd;
          % - now the early etc. peri-event firing rates
          % a row array listing the average (across bins) firing rate of each
          % of the bursts
          tmpFr=mean(hst(pethBinIx{1},:),1)/ap.pethBinW*1000;
          mnPethEarlyFr(1,colIx,masterSliceIx)=nanmean(tmpFr);
          stdPethEarlyFr(1,colIx,masterSliceIx)=nanstd(tmpFr);

          tmpFr=mean(hst(pethBinIx{2},:),1)/ap.pethBinW*1000;
          mnPethIntermedFr(1,colIx,masterSliceIx)=nanmean(tmpFr);
          stdPethIntermedFr(1,colIx,masterSliceIx)=nanstd(tmpFr);

          tmpFr=mean(hst(pethBinIx{3},:),1)/ap.pethBinW*1000;
          mnPethLateFr(1,colIx,masterSliceIx)=nanmean(tmpFr);
          stdPethLateFr(1,colIx,masterSliceIx)=nanstd(tmpFr);
        end
      end
    end % for: curExpChanIx
  end
end

% - compute CVs
cvBurstLen=stdBurstLen./mnBurstLen;
cvSilentPerLen=stdSilentPerLen./mnSilentPerLen;
cvBurstStrength=stdBurstStrength./mnBurstStrength;
cvBurstPeak=stdBurstPeak./mnBurstPeak;
cvPethEarlyFr=stdPethEarlyFr./mnPethEarlyFr;
cvPethIntermedFr=stdPethIntermedFr./mnPethIntermedFr;
cvPethLateFr=stdPethLateFr./mnPethLateFr;

% ----- aftermaths:
% - warn of missing and problematic files
if ~isempty(filesNotFound.fpEtsl)
  msg='The following FP/MP results files (*res*.mat) are listed as useable in the master list but could not be found:';
  warndlg(cat(2,{msg},filesNotFound.fpEtsl));
  warning(msg);
  disp(strvcat(filesNotFound.fpEtsl{:}))
end
if ~isempty(filesEmptyEtsl)
  msg='The following FP/MP results files (*res*.mat) were found but had empty time stamp lists (=no bursts). This is legal, but you may want to check nonethless';
  warndlg(cat(2,{msg},filesEmptyEtsl));
  warning(msg);
  disp(strvcat(filesEmptyEtsl{:}))
end
if ~isempty(filesNotFound.spxTsl)
  msg='the following spikes results files (*SPX_res*.mat) could not be found:';
  warndlg(cat(2,{msg},filesNotFound.spxTsl));
  warning(msg);
  disp(strvcat(filesNotFound.spxTsl{:}))
end
if ~isempty(filesEmptyTsl)
  msg='The following spikes results files (*SPX_res*.mat) were found but had empty time stamp lists. This is legal, but you may want to check nonethless';
  warndlg(cat(2,{msg},filesEmptyTsl));
  warning(msg);
  disp(strvcat(filesEmptyTsl{:}))
end
if ~isempty(filesNotFound.raw)
  msg='The following raw/preconditioned data files (*.mat or *.abf)  could not be found:';
  warndlg(cat(2,{msg},filesNotFound.raw));
  warning(msg);
  disp(strvcat(filesNotFound.raw{:}))
end

% - *** generate normalized versions of all dependent vars ***
for k=1:nMasterDepPar
  y=eval([masterDepPar{k,1} ';']);
  switch masterDepPar{k,2}
    case 1
      % scalar results
      y=y./repmat(y(1,indepParNormIx,:),[1 indepParNLevel 1]);
    case 2
      % peths
      warning('normalization of peth not yet implemented');
    otherwise
      error('normalization of dep var other than scalar or peth not yet implemented');
  end
  eval([masterDepPar{k,1} '_norm=y;']);
end

% save results
% - rename ap and indepPar and save them too, along with other vital information
ap_s=ap;
indepPar_s=indepPar;
tmp=strcat(masterDepPar(:,1),'_norm');
save([ap.resPath '\' ap.resFn],...
  'ap_s','indepPar_s','masterDepPar','expChanName','fileInfo',...
  'drug1_applicRank','undef1','undef2',...
  'Etsl','SilentEtsl','Tsl',...
  'psdF','psDensMn','psDensStd','autoC',...
  'pethBin','pethMn','pethStd',...
  'indepParLevel','indepParNLevel','indepParNormIx',...
  masterDepPar{:,1},tmp{:});

  