function indepPar=tslbatch_readwsdata(fn,indepPar)
% ** function indepPar=tslbatch_readwsdata(fn,indepPar) 
% reads data from an excel worksheet, checks the data and adds data columns
% as field 'd' to struct 'indepPar'

nPar=length(indepPar);
[num,cha,raw] = xlsread(fn,1);
% retrieve the 'header line' (the first element of cha)
header=cha{1};
% headers of the columns
colHead=cha(2,:);
% remove both from raw
raw(1:2,:)=[];
% now the size of raw reflects pure content of the spreadsheet
[n1,n2]=size(raw);

% --- perform some global checks
if n2~=nPar
  errordlg(['the master file contains too few or too many columns (this may be due to a ''stray cell'' inadvertently filled with a value)']);
  error('see error window for description of problem');
end
% if num and raw have different number of lines we probably have a stray
% value in a line
if n1~=size(num,1)
  errordlg(['the master file contains one or more incomplete lines (this may be due to a ''stray cell'' inadvertently filled with a value)']);
  error('see error window for description of problem');
end

for g=1:nPar
  switch indepPar(g).dType
    case 'char'
      if isnumeric(raw(:,g))
        error(['parameter ' indepPar(g).name ' should be char but is numeric']);
      else
        indepPar(g).d=raw(:,g);
      end
    otherwise
      tmp=cat(1,raw{:,g});
      if ~isnumeric(tmp)
        error(['parameter ' indepPar(g).name ' should be numeric but is not']);
      else
        indepPar(g).d=tmp;
      end
  end
end


% ------------ some checks on contents -----------------------------------
parNm={indepPar.name};
% some indices into indepPar needed throughout
expIDIx=strmatch('expID',parNm);
isValidIx=strmatch('isValid',parNm);
fileNameIx=strmatch('fileName',parNm);
chanNameIx=strmatch('chanName',parNm);
drug1_concIx=strmatch('drug1_conc',parNm);
drug1_applicRankIx=strmatch('drug1_applicRank',parNm);

% check whether date part of expID properly reflects file names
tmpStr1=strvcat(indepPar(expIDIx).d);
tmpStr2=strvcat(indepPar(fileNameIx).d);
tmp=strcmp(cellstr(tmpStr1(:,6:15)),tmpStr2(:,1:10));
if ~all(tmp)
  errordlg({'File names'; ' '; tmpStr2(~tmp,:); ' '; 'are associated with mismatching experimental IDs.'});
  error('see error dialog');
end

% drug appl rank must be ints
if sum(double(uint8(indepPar(drug1_applicRankIx).d))-indepPar(drug1_applicRankIx).d)
  error('drug 1 application ranks must be integers 0, 1, 2, etc.');
end














