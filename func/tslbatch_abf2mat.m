function tslbatch_abf2mat(indepPar,ap)
% ** function tslbatch_abf2mat(indepPar,ap)
% reads abf files in list, filters and downsamples them according to
% entries in ap and saves them using abfmerge2mat as *.mat files

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

% names & number of parameters in indepPar
parNm={indepPar.name};
nPar=length(indepPar);

% some indices into indepPar needed throughout
expIDIx=strmatch('expID',parNm);
isValidIx=strmatch('isValid',parNm);
fileNameIx=strmatch('fileName',parNm);
chanNameIx=strmatch('chanName',parNm);

% - check ap & set flags

nDepPar=size(ap.depPar,1);

% -------------------------------------------------------------------------
%                PART I: purges of observations
% -------------------------------------------------------------------------
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
% *** first thing to do here: restructure entries of indepPar(expIDIx).d such
% that they are identical to the directory names
dDir=char(indepPar(expIDIx).d);
% watch out, dDir is a char array
dDir=dDir(:,[6:16 1:4]);
% unique directory names
[uDir,ix1,ix2]=unique(dDir,'rows');
nDir=size(uDir,1);

disp(['** ' int2str(nDir) ' directories'])

% -------------------------------------------------------------------------
%                PART III: collect data
% -------------------------------------------------------------------------
for g=1:nDir
  % find all unique channel names for current dir and convert only these
  uChan=unique(indepPar(chanNameIx).d(ix2==g));
  % find all GAP-FREE *.abf files in current directory and put them on the list
  datDir=[ap.rootPath uDir(g,:)];
  s=dir([datDir '\*.abf']);
  for fIx=1:numel(s)
    [nada,niente,fi]=abfload([datDir '\' s(fIx).name],'info');
    if fi.nOperationMode==3
      abfmerge2mat({s(fIx).name(1:end-4)},'dDir',datDir,'sampFreq',ap.sampleFreq,...
        'cFreq',[ap.hiCFreq ap.loCFreq],'ch',uChan);
    end
  end    
end

