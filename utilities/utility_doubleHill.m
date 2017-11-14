clear

% ******************* USER INPUT ******************************************
% §§§ mat file name (output of tslbatch)
mFn='G:\berthold\Diazepam_SZ\diazepam_ALLE';
% §§§ concentrations IN ASCENDING ORDER
curIpVal=[0 0.005 0.01 0.05 0.25 0.5 0.75 1 6.25 12.5 25 50 100];
% parameter of interest: firing rate, normalized (to control)
par={'eventRate_norm'};
% par={'relTimeInBurst_norm'};
% **** job to be performed & specific options ***

% graphics options
labelscale('fontSz',12,'scaleFac',1,'lineW',2,'markSz',10);
% §§§ name & full path to figure to be saved
figName='G:\berthold\Diazepam_SZ\Dia_doubleHill';

% printas=[];  % don't save graphics
% printas='-djpeg97'; % save as jpeg ('working style')
printas='-dpsc2'; % save as color post script (nice finish)

% ******************* END USER INPUT **************************************

etslconst;
nIpVal=numel(curIpVal);
% load *.mat data file
load(mFn);

% index to columns in all variables
[nix,columnIx]=intersect(indepParLevel,curIpVal);
if ~isequal(nix,curIpVal)
  errordlg({'variable ''curIpVal'' contains values that do not exist in the data.',...
    'Available values: ', num2str(indepParLevel)})
  return
end

d=eval(par{1});
d(~isfinite(d))=nan;
if true
  % §§§ fit data using individual experimental data (need to be reshaped to 2D)
  d=permute(d,[3 2 1]);
  curIpVal=repmat(curIpVal,size(d,1),1);
  d=d(:,columnIx,:);
  curIpVal=curIpVal(:);
  d=d(:);
  goodIx=(isfinite(d));
  ds.crf=[curIpVal(goodIx),d(goodIx)];
else
  % fit average data
  d=nanmean(d,3)';
  goodIx=(isfinite(d));
  ds.crf=[curIpVal(goodIx)',d(goodIx)];
end


% % §§§ TEST DATA
% ds.crf=[...
%   .005 1;...
%   .01  1;...
%   .05  .8;...
%   .2   .82;...
%   .4   .8;...
%   .8   .81;...
%   1    .8;
%   6    .81;...
%   15   .75;...
%   30   .65;...
%   50   .41;...
%   100  .18];
% 



% ---------------------------------------------------------

nDs=1;
isNormalizedData=true;
invertVal=1;
% §§§
% errType='ci';
errType='sem';


for g=1:nDs
  % get rid of nans
  ix=find(any(~isfinite(ds(g).crf),2));
  if ~isempty(ix)
    disp(['data set # ' int2str(g) ': purging ' int2str(length(ix)) ' bad entries']);
    ds(g).crf(ix,:)=[];
  end
  % compute means & std of dependent var
  % - unitary values of independent var (conc)
  ds(g).uIdepVal=unique(ds(g).crf(:,1));
  % - preallocate mn & std
  ds(g).mn=nan*ds(g).uIdepVal;
  ds(g).st=nan*ds(g).uIdepVal;
  for h=1:numel(ds(g).uIdepVal)
    ix=ds(g).crf(:,1)==ds(g).uIdepVal(h);
    ds(g).n(h)=length(find(ix));
    ds(g).mn(h)=mean(ds(g).crf(ix,2));
    ds(g).st(h)=std(ds(g).crf(ix,2));
    ds(g).sem(h)=ds(g).st(h)/sqrt(ds(g).n(h));
    ds(g).ci(h)=ds(g).sem(h)* -tinv(.05/2,ds(g).n(h)-1)
  end
  % if this is normalized data do not fit to control values because they
  % have zero variablility 
  if isNormalizedData
    ix=isalmost(ds(g).crf(:,1),0,eps);
    disp(['omitting ' int2str(numel(find(ix))) ' control entries from fit']);
    ds(g).crf(ix,:)=[];
  end
  % invert?
  if isfinite(invertVal)
    ds(g).crf(:,2)=invertVal-ds(g).crf(:,2);
  end
  % sort according to conc 
  ds(g).crf=sortrows(ds(g).crf,1);
end

% set default colors & markers in case they were not defined
if ~isfield(ds,'pCol')
  ds(1).pCol='k';
  ds(1).pSymb='o';
  if nDs>1
    ds(2).pCol='r';
    ds(2).pSymb='s';
  end
end

% ------- set up fitting parameters ----------------
if 0
  % ------- Hill fit
  ft_ = fittype(['v*conc^n/(k^n + conc^n)'] ,...
    'independent',{'conc'},...
    'coefficients',{'v','n','k'});
  fo_ = fitoptions('method','NonlinearLeastSquares','upper',[1 10000 10000]);
  st_ = [1 1 10];
  nFitPar=3;
else
  ft_ = fittype(['v1*conc^n1/(k1^n1 + conc^n1) + v2*conc^n2/(k2^n2 + conc^n2)'] ,...
    'independent',{'conc'},...
    'coefficients',{'v1','n1','k1','v2','n2','k2'});
  fo_ = fitoptions('method','NonlinearLeastSquares',...
    'upper',[.5    100   .1        1  10 100],...
    'lower',[.01  .5     .01      .3  1    10]);
  st_ = [.2 1 .02     .6 1 30];
  nFitPar=6;
end


% ------- fit ----------------
for g=1:nDs
  set(fo_,'Startpoint',st_);
  [ds(g).fit,ds(g).gof]=fit(ds(g).crf(:,1),ds(g).crf(:,2),ft_ ,fo_);
  ds(g).crf(:,3)=ds(g).fit(ds(g).crf(:,1));
  % extract fit parameters to replace guessed start parameters by the ones 
  % just calculated and to determine degrees of freedom 
  tmpPar=coeffvalues(ds(g).fit);
  ds(g).df=size(ds(g).crf(:,1),1)-numel(tmpPar);  
  ds(g).hillC.vMax=tmpPar(1);
  ds(g).hillC.n=tmpPar(2);
%  ds(g).hillC.k=tmpPar(3); 
end


for g=1:nDs
  if isfinite(invertVal)
    ds(g).crf(:,[2 3])=invertVal-ds(g).crf(:,[2 3]);  
  end
end

figure(1), clf, hold on

for g=1:nDs
  % produce fit over whole range of indep var 
  plotX=ds(g).crf([1 end],1);
  plotX=logspace(log10(ds(g).crf(1,1)),log10(ds(g).crf(end,1)),200);
  plotY=ds(g).fit(plotX);
  % fits just produced must be inverted
  if isfinite(invertVal)
    plotY=invertVal-plotY;
  end
  ph=errorbar(ds(g).uIdepVal,ds(g).mn,ds(g).(errType),ds(g).pSymb);
  set(ph,'color',ds(g).pCol,'markerfacecolor',ds(g).pCol);
  ph=plot(plotX,plotY,'-');
  set(ph,'color',ds(g).pCol);
end
nicexyax(40)
set(gca,'xscale','log')

% - legend 
% collect parameters down along columns (ds along rows) 
% tmpPar=repmat(nan,nFitPar+1,nDs);
% for g=1:nDs
%   % §§ doesn't work with linear - adjust to number of pars fitted
%   tmpPar(1:3,g)=coeffvalues(ds(g).fit)';
%   % append Radj
%   tmpPar(4,g)=ds(g).gof.adjrsquare;
% end
% 
% txt={...
%   ['EC_{50}: ' num2str(tmpPar(3,:),'%3.2f ') ' nM'],...
%   ['Limit: ' num2str(tmpPar(1,:),'%3.2f ')],...
%   ['R_{adj}: ' num2str(tmpPar(4,:),'%3.2f ')],...
%   };
% % dummy handle
% dummHand=line(nan,nan,'linestyle','none');
% wp.legH=legend(repmat(dummHand,numel(txt),1),txt);
% nicexyax



% ----- printing -----
if ~isempty(printas), 
  print(printas,figName,'-r400'); 
end