function indepPar=tslbatch_defdata(ap)
% ** function function indepPar=tslbatch_defdata(ap)
% generates structure array dsPar with as many elements as columns of a
% worksheet (*.xls) to be read by read_wsdata.m. The definition
% of the columns of the worksheet is in the body of the function. The input
% arg allow specification of independent and dependent variables as well
% as values of certain parameters which shall function as exclusion criteria

% *************************************************************************
% Below is the definition of columns in the worksheet. Don't change.
% *************************************************************************
g=1;
indepPar(g).name='expID';
indepPar(g).dType='char';

g=g+1;
indepPar(g).name='fileName';
indepPar(g).dType='char';

g=g+1;
indepPar(g).name='chanName';
indepPar(g).dType='char';

g=g+1;
indepPar(g).name='drug1_conc';
indepPar(g).dType='double';

g=g+1;
indepPar(g).name='drug1_applicRank';
indepPar(g).dType='double';

g=g+1;
indepPar(g).name='drug2_conc';
indepPar(g).dType='double';

g=g+1;
indepPar(g).name='isValid';
indepPar(g).dType='double';

g=g+1;
indepPar(g).name='ageAtPrep';
indepPar(g).dType='double';

g=g+1;
indepPar(g).name='prepID';
indepPar(g).dType='double';

g=g+1;
indepPar(g).name='daysInV';
indepPar(g).dType='double';

g=g+1;
indepPar(g).name='typeCulture';
indepPar(g).dType='char';

g=g+1;
indepPar(g).name='genotype';
indepPar(g).dType='char';

g=g+1;
indepPar(g).name='undef1';
indepPar(g).dType='double';

g=g+1;
indepPar(g).name='undef2';
indepPar(g).dType='double';

g=g+1;
indepPar(g).name='comment';
indepPar(g).dType='char';


% -------- complement indepPar and check ds ------------------------
nPar=length(indepPar);
charIx=strmatch('char',strvcat(indepPar.dType));
numIx=setdiff(1:nPar,charIx);
% append fields to indepPar and set values based on ap.indepPar 
[indepPar.use]=deal(0);
% set .normVal to nan or '' depending on type
[indepPar(charIx).normVal]=deal('');
[indepPar(numIx).normVal]=deal(nan);
ix=ismember(ap.indepPar(:,1),{indepPar.name});
if isempty(ix)
  errordlg({'independent parameter(s) is/are misspelled or nonexistent'});
  error('see error window for description of problem');
elseif ~all(ix)
  errordlg({'independent parameter(s)', ' ', ap.indepPar{~ix,1}, ' ','is/are misspelled or nonexistent'});
  error('see error window for description of problem');  
elseif numel(ix)>2
  errordlg(['too many independent analysis parameters specified (presently, max. of two)']);
  error('see error window for description of problem');
else
  [nix,ix1,ix2]=intersect({indepPar.name},ap.indepPar(:,1));
  [indepPar(ix1).use]=deal(1);
  [indepPar(ix1).normVal]=deal(ap.indepPar{ix2,2});
end

% same story for ap.subsetPar
[indepPar(numIx).min]=deal(nan);
[indepPar(numIx).max]=deal(nan);
[indepPar(numIx).purgeVal]=deal(nan);
[indepPar(charIx).min]=deal('');
[indepPar(charIx).max]=deal('');
[indepPar(charIx).purgeVal]=deal('');
ix=ismember(ap.subsetPar(:,1),{indepPar.name});
if isempty(ix)
  errordlg({'subset parameter(s) is/are misspelled or nonexistent'});
  error('see error window for description of problem');
elseif ~all(ix)
  errordlg({'subset parameter(s)', ' ', ap.subsetPar{~ix,1}, ' ','is/are misspelled or nonexistent'});
  error('see error window for description of problem');  
else
  [nix,ix1,ix2]=intersect({indepPar.name},ap.subsetPar(:,1));
  [indepPar(ix1).min]=deal(ap.subsetPar{ix2,2});
  [indepPar(ix1).max]=deal(ap.subsetPar{ix2,3});  
  [indepPar(ix1).purgeVal]=deal(ap.subsetPar{ix2,4});  
end



