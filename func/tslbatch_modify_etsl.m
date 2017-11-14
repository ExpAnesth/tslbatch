function tslbatch_modify_etsl(indepPar,ap)
% ** function tslbatch_modify_etsl(indepPar,ap)
% loops through etsl and shifts time zero of events 

% UNDER CONSTRUCTION!!

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

% interval in which to determine tZero
wp.tZeroFindIntv=[-5 500];

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
recSiteIx=strmatch('recSite',parNm);
undef1Ix=strmatch('undef1',parNm);
undef2Ix=strmatch('undef2',parNm);

etslconst;

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
    warndlg({'In ';' '; strvcat(indepPar(expIDIx).d{tmpIx}); ' '; 'at least one control recording (defined as [drug1]=0) has an application rank different from 0 - make sure that this is correct'});
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

% important: have a 1D cell array with as many elements as the results
% variables above have 'slices': it contains strings identifying the
% nUExpChan 'experiments'
expChanName=cell(1,nUExpChan);
% we also need rudimentary file info
fileInfo.dDir='';
fileInfo.fileName='';
fileInfo.chName='';
fileInfo.recSite=nan;
fileInfo.recTime=nan;
fileInfo.si=nan;
fileInfo.winBuCutout=ap.winBuCutout;
fileInfo=repmat(fileInfo,[1 indepParNLevel nUExpChan]);

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
      % ** set flags for the different analyses

      % file name without extension as listed in master table
      curFnBase=indepPar(fileNameIx).d{curExpChanIx(i)};
      % name of FP results file
      curFnFpRes=[curFnBase '_' curDeblankChName '_FP_res'];
      % name of MP results file
      curFnMpRes=[curFnBase '_' curDeblankChName '_MP_res'];
      disp(['*** ' curFnBase '_' curDeblankChName '...']);

      % collect file name without extension & other info no matter whether
      % the data were found or not
      fileInfo(1,colIx,masterSliceIx).fileName=curFnBase;
      fileInfo(1,colIx,masterSliceIx).dDir=datDir;
      fileInfo(1,colIx,masterSliceIx).chName=curDeblankChName;   
      fileInfo(1,colIx,masterSliceIx).recSite=indepPar(recSiteIx).d(curExpChanIx(i));
      % also, collect all independent parameters to be saved in a format
      % identical to the scalar results variables
      drug1_applicRank(1,colIx,masterSliceIx)=indepPar(drug1_applicRankIx).d(curExpChanIx(i));
      recSite(1,colIx,masterSliceIx)=indepPar(recSiteIx).d(curExpChanIx(i));
      undef1(1,colIx,masterSliceIx)=indepPar(undef1Ix).d(curExpChanIx(i));
      undef2(1,colIx,masterSliceIx)=indepPar(undef2Ix).d(curExpChanIx(i));      
      % -------------------------------------------------------------------
      % ---------------------- reading data -------------------------------
      % -------------------------------------------------------------------
      
      % ------- load general information ------------      
      % essential information on the recording can be extracted from any of
      % the files pertaining to a recording, but as the odd file may be
      % missing try to retrieve this information here, independent of the
      % data proper.
      fList={[datDir '\' curFnFpRes '.mat'],...
        [datDir '\' curFnMpRes '.mat']};
      fExist=[exist(fList{1},'file'),exist(fList{2},'file')];
      % pick FP/MP file for header info
      if any(fExist)
        load(fList{find(fExist,1)},'head');
        % recording time in s
        fileInfo(1,colIx,masterSliceIx).recTime=head.ds.fileInfo.recTime;
        % sampling interval of downsampled (interpolated) data used here
        fileInfo(1,colIx,masterSliceIx).si=1e6/ap.sampleFreq;
      end
      % generate empty vars
      etsl=[];
      silentEtsl=[];
      tsl=[];

      % ------- load FP/MP burst results ------------
      if (isValid==wp.useFP || isValid==wp.useFPSPX)
        fExist=fExist(1:2);
        if any(fExist)
          if fExist(1)
            isExtraRec=true;
            load([datDir '\' curFnFpRes '.mat']);
          else
            load([datDir '\' curFnMpRes '.mat']);
            isExtraRec=false;
          end
          % detach etsl, silentEtsl and stats
          etsl=bu.etsl;
          silentEtsl=bu.silentEtsl;
          if isempty(etsl)
            warning('no FP or MP bursts detected');
            filesEmptyEtsl{end+1}=curFnBase;
          end
        else
          % §§ there's an ambiguity here: it cannot be resolved whether the
          % file looked for but not found is a FP or a MP one. Collect them
          % in the list of missing FP files for now.
          filesNotFound.fpEtsl{end+1}=curFnFpRes;
          etsl=[];
          silentEtsl=[];
        end
      end
      
      % -------------- raw or preconditioned FP data ----------------------
      % no need to check for isValid, has been done above
      
      isExistentFile=[exist([datDir '\' curFnBase '.mat'],'file'),...
        exist([datDir '\' curFnBase '.abf'],'file')];
      if ~any(isExistentFile)
        filesNotFound.raw{end+1}=curFnBase;
        doBIA=false;
        doAC=false;
        doSA=false;
        doSTA=false;
        % no need to do anything to collected results variables because
        % they are preallocated with nans (all other)
      else
        % ** file loading logics:
        % - if mat raw file exists load it
        % - if only abf file exists load it and preprocess data (time-
        % consuming)
        if isExistentFile(1)
          [d,si,fi]=matDload([datDir '\' curFnBase '.mat'],'channels',{curDeblankChName});
          % clear si because it is not meant to be used (use fi.si)
          clear si
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
          maxNExcPts=diff(cont2discrete(wp.tZeroFindIntv,1000/ap.sampleFreq,'intv',1))+1;
          % - lengths of bursts in ticks
          buL=cont2discrete(etsl(:,etslc.durCol),1000/ap.sampleFreq,'intv',1);
          % - beginnings of bursts in ticks
          buB=cont2discrete(etsl(:,etslc.tsCol),1000/ap.sampleFreq,'intv',1);
          % - cutout times according to wp.tZeroFindIntv in ticks
          buC=cont2discrete(repmat(etsl(:,etslc.tsCol),1,2)+...
            repmat(wp.tZeroFindIntv,nBu,1),1000/ap.sampleFreq,'intv',1);
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
          sum(etsl(:,[etslc.tsCol etslc.durCol]),2)+ap.winBuCutout(2)< diff(head.ds.fileInfo.recTime)*1000;
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
        % ********* below, commented: choice of end points as of March 14,
        % 2011, namely, up to the end of the bursts
        % buE=cont2discrete(sum(etsl(:,[etslc.tsCol etslc.durCol]),2)+ap.winBuCutout(2),1000/ap.sampleFreq,'intv',1);
        % current version: ap.ppSeg
        buE=cont2discrete(etsl(:,etslc.tsCol)+ap.winBuCutout(1)+ap.ppSeg,1000/ap.sampleFreq,'intv',1);
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

      
    end % for: curExpChanIx
  end
end


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
if ~isempty(filesNotFound.raw)
  msg='The following raw/preconditioned data files (*.mat or *.abf)  could not be found:';
  warndlg(cat(2,{msg},filesNotFound.raw));
  warning(msg);
  disp(strvcat(filesNotFound.raw{:}))
end

  