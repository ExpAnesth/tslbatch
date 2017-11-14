function tslbatch_procdata(indepPar,ap)
% ** function tslbatch_procdata(indepPar,ap)
% core processing function of tslbatch workhorse. See the documentation and
% end of this m-file.

% -------------------------------------------------------------------------
%             PART 0: preliminaries: auxil pars & checks of ap
% -------------------------------------------------------------------------
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

etslconst;

% - check ap & set flags
% if this flag is true the original value shall be used
doUseOrig_minSilentPerDur=false;
if ischar(ap.minSilentPerDur)
  if strcmp(ap.minSilentPerDur,'original')
    doUseOrig_minSilentPerDur=true;
  else
    error('''ap.minSilentPerDur'' must be a number or string ''original''');
  end
end
doUseOrig_minSilentPerNEv=false;
if ischar(ap.minSilentPerNEv)
  if strcmp(ap.minSilentPerNEv,'original')
    doUseOrig_minSilentPerNEv=true;
  else
    error('''ap.minSilentPerNEv'' must be a number or string ''original''');
  end
end
% doUseOrig_maxNormalBurstLength=false;
% % this is a tribute to the fact that this parameter has been introduced
% % quite late
% if ~isfield(ap,'maxNormalBurstLength')
%   ap.maxNormalBurstLength=inf;
% end
% if ischar(ap.maxNormalBurstLength)
%   if strcmp(ap.maxNormalBurstLength,'original')
%     doUseOrig_maxNormalBurstLength=true;
%   else
%     error('''ap.maxNormalBurstLength'' must be a number or string ''original''');
%   end
% end


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


% Make sure that all dependent parameters requested are part of the
% staple computed here
% The second column of masterDepPar indicates whether the dependent
% variable is scalar (1), or sth else (currently not used); the third
% column indicates whether the parameter is a FP-related (1) or a
% SPX-related one (2) - this information will be used later for the
% generation of summary plots
masterDepPar={...
  'eventRate'            , 1, 2;... % firing rates
  'cvISI'                , 1, 2;...% cv of inter-spike-interval
  'fractionEvInBurst'    , 1, 2;... % fraction of spx in bursts
  'relTimeInBurst'       , 1, 2;... % relative time spent in bursts ('active' time)
  'mnNEvPerBurst'        , 1, 2;... % mean of number of events in burst
  'stdNEvPerBurst'       , 1, 2;... % std of number of events in burst
  'cvNEvPerBurst'        , 1, 2;... % cv of number of events in burst
  'mnIntraBurstEvRate'   , 1, 2;... % mean of intra-burst event rate
  'stdIntraBurstEvRate'  , 1, 2;... % std of intra-burst event rate
  'cvIntraBurstEvRate'   , 1, 2;... % cv of intra-burst event rate
  'mnBurstLen'           , 1, 2;... % mean length of bursts
  'mdBurstLen',           1, 2;... % median of length of bursts
  'stdBurstLen'          , 1, 2;... % std of length of bursts
  'cvBurstLen'           , 1, 2;... % cv of length of bursts
  'asBurstLen',           1, 2;... % asymmetry of length of bursts: (95th percentile-median)/(median-5th percentile)
  'burstRate'            , 1, 2;... % burst rate
  'mnSilentPerLen'       , 1, 2;... % mean length of silent period
  'mdSilentPerLen',       1, 2;... % median of length of silent period
  'stdSilentPerLen'      , 1, 2;... % std of length of silent period
  'cvSilentPerLen'       , 1, 2;... % cv of length of silent period
  'asSilentPerLen',       1, 2;... % asymmetry of length of silent period
  'mnPethEarlyFr'        , 1, 2;... % mean of early firing rate
  'stdPethEarlyFr'       , 1, 2;... % std of early firing rate
  'cvPethEarlyFr'        , 1, 2;... % cv of early firing rate
  'mnPethIntermedFr'     , 1, 2;... % mean of intermediate firing rate
  'stdPethIntermedFr'    , 1, 2;... % std of intermediate firing rate  
  'cvPethIntermedFr'     , 1, 2;... % cv of intermediate firing rate
  'mnPethLateFr'         , 1, 2;... % mean of late firing rate  
  'stdPethLateFr'        , 1, 2;... % std of late firing rate
  'cvPethLateFr'         , 1, 2;... % cv of late firing rate    
};

tmp=setdiff(ap.depPar(:,1),masterDepPar(:,1));
if ~isempty(tmp)
  errordlg({'Analysis parameters'; ' '; strvcat(tmp); ' '; 'do not exist (check spelling - NO BLANKS ARE ALLOWED).'; ' '; 'Here is the full list of legal parameters:'; ' '; strvcat(masterDepPar(:,1))});
  error('see error dialog');
end
nMasterDepPar=size(masterDepPar,1);
nDepPar=size(ap.depPar,1);

% -------------------------------------------------------------------------
%                PART I: purges of observations
% -------------------------------------------------------------------------
% §§ convert text parameters like cultureType that may eventually be used
% as indep par to numbers


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


% ---------- preallocation of processed results structure --------------

% --- dimensions of cell array containing RAW values (tsl):
%   <one row> | indep par value | <unique combin experiment, channel>
rawTemplate=cell([1 indepParNLevel nUExpChan]);
% note UPPERCASE
Tsl=rawTemplate;
Etsl=rawTemplate;
SilentEtsl=rawTemplate;

% --- dimensions of histogram variables:
%   peri-event time bins| indep par value | <unique combin experiment, channel>
pethTemplate=repmat(hstTemplate,[1 indepParNLevel nUExpChan]);
pethMn=pethTemplate;
pethStd=pethTemplate;

% --- dimensions of SCALAR values like firing rates:
%   <one row> | indep par value | <unique combin experiment, channel>
scalarTemplate=repmat(nan,[1 indepParNLevel nUExpChan]);


% - firing rates
eventRate=scalarTemplate;
% - cv of ISI
cvISI=scalarTemplate;
% - fraction of spx in bursts
fractionEvInBurst=scalarTemplate;
% - relative time spent in bursts ('active' time)
relTimeInBurst=scalarTemplate;
% - mean, std & cv of number of events in burst
mnNEvPerBurst=scalarTemplate;
stdNEvPerBurst=scalarTemplate;
cvNEvPerBurst=scalarTemplate;
% - mean, std & cv of intra-burst event rate
mnIntraBurstEvRate=scalarTemplate;
stdIntraBurstEvRate=scalarTemplate;
cvIntraBurstEvRate=scalarTemplate;
% - mean, std, cv, etc. of length of bursts
mnBurstLen=scalarTemplate;
mdBurstLen=scalarTemplate;
stdBurstLen=scalarTemplate;
cvBurstLen=scalarTemplate;
asBurstLen=scalarTemplate;
% - burst rate
burstRate=scalarTemplate;
% - mean, std, cv, etc. of length of silent period
mnSilentPerLen=scalarTemplate;
mdSilentPerLen=scalarTemplate;
stdSilentPerLen=scalarTemplate;
cvSilentPerLen=scalarTemplate;
asSilentPerLen=scalarTemplate;
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
% also export application rank order (may be useful)
drug1_applicRank=scalarTemplate;

% important: have a 1D cell array with as many elements as the results
% variables above have 'slices': it contains strings identifying the
% nUExpChan 'experiments'
expChanName=cell(1,nUExpChan);

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
% container for files not found
filesNotFound=cell(0);
% container for files without prior burst detection
filesNoBurstAnalysis=cell(0);
% container for files with empty tsl (=zero firing rate)
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
    masterSliceIx=masterSliceIx+1;
    % fill expChanName: experiment ID + file name
    expChanName{masterSliceIx}=[indepPar(expIDIx).d{curExpChanIx(1)},...
      ', ', indepPar(chanNameIx).d{curExpChanIx(1)}];
    % --- for given channel, loop over files (=levels of indep par)
    for i=1:numel(curExpChanIx)
      % value of level
      curLevel=indepPar(aIndepIx).d(curExpChanIx(i));
      % index to column
      colIx=curLevel==indepParLevel;
      % by default, of course do burst analysis; however, depending on
      % data, change that
      doBurstAnalysis=true;
      % ---------------- read atsl/(e)tsl ------------------------------
      % deblanked channel name
      dbch=indepPar(chanNameIx).d{curExpChanIx(i)};
      dbch=dbch(~isspace(dbch));
      tmpFnBody=[indepPar(fileNameIx).d{curExpChanIx(i)} '_' dbch '_spx_res'];
      isExistentFile=[exist([datDir '\' tmpFnBody '.mat.txt'],'file'),...
        exist([datDir '\' tmpFnBody '.mat'],'file')];
%       % find out whether ap_bu exists in matfile: if that is not so either
%       % burst detection has not been performed or it has been done with an
%       % earlier version of burstdetgui. In either case, don't use the
%       % matfile.
%       if isExistentFile(2)==2
%         fileContent=whos('-file',[datDir '\' tmpFnBody '.mat']);
%         if isempty(strmatch('ap_bu',{fileContent.name}))
%           isExistentFile(2)=0;
%         end
%       end
      % ** file usage rules:
      % - if both files are present use txt file
      % - if only either file present use it
      % - if neither file found, collect name and spit out in warning message at end
      % - if only mat file is present and it does not contain burst
      % parameters issue a warning at the end
      if ~any(isExistentFile)
        filesNotFound{end+1}=tmpFnBody;
        % no need to do anything to collected results variables because
        % they are preallocated with empty matrix (cell array of tsls) or 
        % nans (all other)
      else
        if isequal(isExistentFile,[2 0]) || isequal(isExistentFile,[2 2])
          % if a *.txt atsl file exists it is assumed that burst analysis
          % had been performed on the data so there is no need to check for
          % it
          [head,atsl,etsl,tsl]=read_atsl([datDir '\' tmpFnBody '.mat.txt']);
          % ** create silent etsl: realize that because the txt file contains only 
          % active periods (UP states) silent periods cannot be recreated
          % in the same way as from raw data
          % beginning of burst (t+1) - end of burst (t)
          silentEtsl=etsl(2:end,etslc.tsCol)-sum(etsl(1:end-1,[etslc.tsCol etslc.durCol]),2);
        else
          load([datDir '\' tmpFnBody '.mat']);
          % detach tsl, etsl and silentEtsl
          etsl=bu.etsl;
          silentEtsl=bu.silentEtsl;
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
          clear bu evt
          % generate generic 'head' using define_atsl.m
          [atslDef,tmpHead]=define_atsl;
          % assign common values needed below
          tmpHead.rectime=diff(head.ds.fileInfo.recTime);
          % *** now check whether etsl is empty - if that is so no burst
          % analysis had been performed, which means that none should be
          % done here, either; furthermore, most burst-related analysis
          % parameters (ap_bu.whatever) will not exist
          if ~isempty(etsl)
            try
              tmpHead.burst_interval1=ap_bu.maxIEI_init;
              tmpHead.burst_interval2=ap_bu.maxIEI_tail;
              tmpHead.burst_minev=ap_bu.minNEvPerBurst;
              tmpHead.silentper_dur=ap_bu.minPreBurstSilentPerLen;
              tmpHead.silentper_maxev=ap_bu.maxPreBurstSilentPerNEv;
            catch
              warning('ap_bu not present in mat file (earlier version of burstdetgui');
            end
          end
          origHead=head;
          head=tmpHead;
          clear tmpHead ap_bu;
        end
        % collect tsl and determine how to proceed:
        if isempty(tsl)
          % we don't want to insert an empty matrix where there is already an
          % array of NaNs, so leave Tsl alone
          warning('current tsl is empty');
          filesEmptyTsl{end+1}=tmpFnBody;
          doBurstAnalysis=false;
        else
          Tsl{1,colIx,masterSliceIx}=tsl;
        end
        % if etsl is empty do not perform burst analysis
        if isempty(etsl)
          filesNoBurstAnalysis{end+1}=tmpFnBody;
          doBurstAnalysis=false;
        end

%         % §§§****&&&$$§§§§§******§§§****&&&$$§§§§§******§§§****&&&$$§§§§§******
%         % pssst! secret, provisional code to produce cutouts from
%         % raw data and tsl...
% 
%         if ~isempty(tsl)
%           % adjust directory and channel names
%           ddi=strrep(origHead.ds.dataPath,'E:\ll\rawData\','D:\_data\STN\');
%           fn=[ddi origHead.ds.dataFn];
%           if ~isempty(strfind(origHead.ap.resFn,'_IN0_'))
%             chacha={'IN 0'};
%           else
%             chacha={'IN 1'};
%           end
%           % pick first 100 spx
%           nPickEvt=min(100,numel(tsl)-1);
%           % interval for cutouts (ms)
%           winEvCutout=[-1 2];
%           % load raw data
%           [d,osi]=abfload(fn,'channels',chacha,'stop',tsl(nPickEvt+1)/1000);
%           % filter
%           d=hifi(d,osi,origHead.ap.hiCFreq);
%           % convert event list of detected events to ticks for cutting out excerpts
%           tmpTsl=cont2discrete(tsl,osi/1e3,'intv',1);
%           tmpwinEvCutout=cont2discrete(winEvCutout,osi/1e3,'intv',1);
%           % cut out
%           evCutout=tsl2exc(d,'idx',{tmpTsl(1:nPickEvt)},'win',tmpwinEvCutout);
%           % target sampling interval in us 
%           siTarget=10;
%           % check whether these are more or less stable waveforms
%           r=evdeal(evCutout,'idx','minmaxpeak');
%           if std(r.maxPeak)/mean(r.maxPeak)<.3
%             % new discrete time axis
%             upsampT=linspace(tmpwinEvCutout(1),tmpwinEvCutout(2),(diff(tmpwinEvCutout)+1)*osi/siTarget);
%             [nix,zeroIdx]=min(abs(upsampT));
%             % upsample cutouts
%             upsampEvCutout=interp1(tmpwinEvCutout(1):tmpwinEvCutout(2),evCutout,upsampT,'spline');
%             % average
%             mnEvCutout=mean(upsampEvCutout,2);
%             % determine max peak again
%             r=evdeal(mnEvCutout,'idx','minmaxpeak');
%             if mean(r.maxPeak)> -mean(r.minPeak)
%               tshift=r.maxPeakT-zeroIdx;
%               nFac=r.maxPeak;
%             else
%               tshift=r.minPeakT-zeroIdx;              
%               nFac=-r.minPeak;
%             end
%             % shift circularly...
%             mnEvCutout=circshift(mnEvCutout,-tshift);
%             % normalize to min/max peak
%             mnEvCutout=mnEvCutout/nFac;
%             % preallocate main results variable if it does not yet exist
%             if ~exist('evExc','var')
%               evExc=repmat(nan*mnEvCutout,[1 indepParNLevel nUExpChan]);
%             end
%             evExc(:,colIx,masterSliceIx)=mnEvCutout;
%             % plot(upsampT,mnEvCutout)
%             % drawnow
%           end
%           % to shorten things...
%           doBurstAnalysis=false;
%         end

%  ***** THIS WILL HAVE TO REPLACE SAVE COMMAND BELOW *****
% save([ap.resPath '\' ap.resFn],'ap_s','indepPar_s','expChanName',...
%   'Tsl','Etsl','SilentEtsl',...
%   'pethBin','pethMn','pethStd',...
%   'indepParLevel','indepParNLevel','indepParNormIx',...
%   'drug1_applicRank',...
%   'evExc',...
%   masterDepPar{:,1},tmp{:});

%         % §§§****&&&$$§§§§§******§§§****&&&$$§§§§§******§§§****&&&$$§§§§§******        

        % ---------------- compute parameters ---------------------
        % - event (spx) rate (must be computed here because a firing rate
        % of zero (=empty tsl) is perfectly valid)
        eventRate(1,colIx,masterSliceIx)=numel(tsl)/head.rectime;
        isi=diff(tsl);
        cvISI(1,colIx,masterSliceIx)=std(isi)/mean(isi);
        drug1_applicRank(1,colIx,masterSliceIx)=indepPar(drug1_applicRankIx).d(curExpChanIx(i));
        
        % if this is a recording 
        % - with zero events 
        % - with empty etsl (i.e. no prior burst analysis had been done)
        % skip all of the rest. etslburstf and tslpeth can deal with empty
        % time stamp lists, but averaging over bursts will fail and the
        % whole shebang is futile anyways. All results variables are
        % preallocated with nans, so there will be no funny gaps.
        if doBurstAnalysis
          % *** impose restraining parameters?
          if doUseOrig_minSilentPerDur
            ap.minSilentPerDur=head.silentper_dur;
          end
          if doUseOrig_minSilentPerNEv
            ap.minSilentPerNEv=head.silentper_maxev;
          end
          [etsl,tmpAtsl,silentEtsl,stats]=etslburstf(...
            tsl,head.burst_interval2,...
            'recLen', head.rectime*1000,...
            'maxIEI_init',head.burst_interval1,...
            'minNEvPerBurst',head.burst_minev,...
            'minSilentPerDur',ap.minSilentPerDur,...
            'maxSilentPerNEv',ap.minSilentPerNEv);
          % copy freshly computed (silent) etsls
          Etsl{1,colIx,masterSliceIx}=etsl;
          SilentEtsl{1,colIx,masterSliceIx}=silentEtsl;
          nEv=numel(tsl);
          nBu=size(etsl,1);

          % - fraction of spx in bursts
          fractionEvInBurst(1,colIx,masterSliceIx)=stats.fractionEvInBurst;
          % - relative time spent in bursts ('active' time)
          relTimeInBurst(1,colIx,masterSliceIx)=stats.relTimeInBurst;
          % -  mean, std & cv of number of events in burst
          mnNEvPerBurst(1,colIx,masterSliceIx)=stats.mnNEvPerBurst;
          stdNEvPerBurst(1,colIx,masterSliceIx)=stats.stdNEvPerBurst;
          % -  mean, std & cv of intra-burst event rate
          mnIntraBurstEvRate(1,colIx,masterSliceIx)=stats.mnIntraBurstEvRate;
          stdIntraBurstEvRate(1,colIx,masterSliceIx)=stats.stdIntraBurstEvRate;
          % -  mean, std & cv length of bursts
          mnBurstLen(1,colIx,masterSliceIx)=stats.mnBurstLen;
          stdBurstLen(1,colIx,masterSliceIx)=stats.stdBurstLen;
          % - burst rate
          burstRate(1,colIx,masterSliceIx)=stats.burstRate;
          % -  mean, std & cv length of silent period
          mnSilentPerLen(1,colIx,masterSliceIx)=stats.mnSilentPerLen;
          stdSilentPerLen(1,colIx,masterSliceIx)=stats.stdSilentPerLen;
          % ** add stats not computed by etslburstf:
          if ~isempty(etsl)
            tmpx=prctile(etsl(:,etslc.durCol),[5 50 95]);
            mdBurstLen(1,colIx,masterSliceIx)=tmpx(2);
            asBurstLen(1,colIx,masterSliceIx)=(tmpx(3)-2*tmpx(2)+tmpx(1))/(tmpx(3)-tmpx(1));
          end
          if ~isempty(silentEtsl)
            tmpx=prctile(silentEtsl(:,etslc.durCol),[5 50 95]);
            mdSilentPerLen(1,colIx,masterSliceIx)=tmpx(2);
            asSilentPerLen(1,colIx,masterSliceIx)=(tmpx(3)-2*tmpx(2)+tmpx(1))/(tmpx(3)-tmpx(1));
          end
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
    end
  end
end

% - compute CVs
cvNEvPerBurst=stdNEvPerBurst./mnNEvPerBurst;
cvIntraBurstEvRate=stdIntraBurstEvRate./mnIntraBurstEvRate;
cvBurstLen=stdBurstLen./mnBurstLen;
cvSilentPerLen=stdSilentPerLen./mnSilentPerLen;
cvPethEarlyFr=stdPethEarlyFr./mnPethEarlyFr;
cvPethIntermedFr=stdPethIntermedFr./mnPethIntermedFr;
cvPethLateFr=stdPethLateFr./mnPethLateFr;

% ----- aftermaths:
% - warn of missing and problematic files
if ~isempty(filesNotFound)
  msg='the following files are listed as useable in the master list but could not be found:';
  warndlg(cat(2,{msg},filesNotFound));
  warning(msg);
  disp(strvcat(filesNotFound{:}))
end
if ~isempty(filesNoBurstAnalysis)
  msg='the following files are listed as useable in the master list but did not contain burst analysis parameters:';
  warndlg(cat(2,{msg},filesNoBurstAnalysis));
  warning(msg);
  disp(strvcat(filesNoBurstAnalysis{:}))
end
if ~isempty(filesEmptyTsl)
  msg='the following files are listed as useable in the master list but had empty time stamp lists (=zero firing rate). This is legal, but you may want to check nonethless';
  warndlg(cat(2,{msg},filesEmptyTsl));
  warning(msg);
  disp(strvcat(filesEmptyTsl{:}))
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
save([ap.resPath '\' ap.resFn],'ap_s','indepPar_s','expChanName',...
  'Tsl','Etsl','SilentEtsl',...
  'pethBin','pethMn','pethStd',...
  'indepParLevel','indepParNLevel','indepParNormIx','masterDepPar',...
  'drug1_applicRank',...
  masterDepPar{:,1},tmp{:});

  
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
% INTERNAL NOTE
% In this version of tslbatch all dependent variables will always be
% computed and saved in their entirety. Furthermore, normalized variants 
% of all variables will be generated by default and saved along the
% non-normalized ones. The user has control only over the variables to be
% plotted.
% To implement a version allowing the user full control over the
% computation do the following:
% - use second column of masterDepPar to differentiate between scalar and
% truly 3D variables like histograms. More precisely, collect all results
% in a cell array with as many elements as dep pars and preallocate
% accordingly. This allows for auotmatization of results collection:
%
%       % place all results in a cell array (in the order defined in ap.depPar) and
%       % assign to separate variables only at end
%       r=cell(nDepPar,1);
%       for g=1:nDepPar
%         switch ap.depPar{g,2}
%           case 1
%             r{g}=scalarTemplate;
%           case 2
%             r{g}=pethTemplate;
%           otherwise
%             error('ask the programmer');
%         end
%       end
%
% - run computations conditionally
