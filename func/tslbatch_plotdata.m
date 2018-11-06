function tslbatch_plotdata(fn,varargin)
% ** function tslbatch_plotdata(fn,varargin)
% plotting routine of tslbatch workhorse. See the documentation.

% -------------------------------------------------------------------------
% Version 2.1.3, November 2018
% (C) Harald Hentschke (University Hospital of Tuebingen)
% -------------------------------------------------------------------------

etslconst;
load(fn)
% if ap is given as input arg use it. Be aware that only a few fields have 
% an effect
if nargin<=1
  ap=ap_s;
else
  ap=varargin{1};
end
% rename original indepPar 
indepPar=indepPar_s;

% re-compute some vars
% *** all SPX scalar parameters
[~,dpIx]=union(ap.depPar(:,1),masterDepPar(cat(1,masterDepPar{:,3})==1,1));
ap.depPar=ap.depPar(sort(dpIx),:);
nDepPar=size(ap.depPar,1);
nExpChanName=numel(expChanName);
indepParNLevel=numel(indepParLevel);

% find out whether wash values are included
washIx=find(indepParLevel==-1);
if ~isempty(washIx) 
  if numel(washIx)>1
    error('only one wash recording per experiment allowed');
  else
    if washIx~=1
      error('wash recordings MUST be assigned a drug1 conc value of -1');
    end
  end
  doWash=true;
  % index to dependent vars excluding wash
  nwIx=washIx+1:indepParNLevel;
else
  doWash=false;
  nwIx=1:indepParNLevel;  
end

% internal settings: 
% - type of error bar
ebType='ci'; % sd
% bin width of inst firing rate
ifrBinW=.25;


% *************************************************************************
% *                    summary figure                                     *
% *************************************************************************
indepParLevel_nw=indepParLevel(nwIx);
% --- layout of figure:
% - treat raw peths independent of scalar results
% - reserve first row of figure for peths
% - semi-explicit assignment of rows and columns
if nDepPar<2
  nCol=2;
  nRow=2;
elseif nDepPar<=6
  nCol=2;
  nRow=ceil(nDepPar/nCol)+1;
else
  nCol=3;
  nRow=ceil(nDepPar/nCol)+1;
end
% the subplot offset for scalar plots
scalarSpOffs=nCol;

fh=figure(1); orient tall
clf
ftag='summary plot';
labelscale('fontSz',8,'scaleFac',1.0,'lineW',1.2,'markSz',6);
tmpScrSz=get(0,'Screensize');
tmpScrSz([1 2])=tmpScrSz([1 2])+tmpScrSz([3 4])*.02;
tmpScrSz([3 4])=tmpScrSz([3 4]).*[.5 .9];
set(fh,'position',tmpScrSz,'tag',ftag,'name',ftag,'color',[.9 .9 1],'numbertitle','off');
cm=coma('jet','n',nExpChanName);

% define simple callbacks for mouse click on idividual data points
callB=cell(nExpChanName,1);
for g=1:nExpChanName
  % callB{g}=['title(''' expChanName{g} ''',''interpreter'',''none'');'];
  callB{g}=['disp(''' expChanName{g} ''');'];
end

% plot peths
pethMnMn=nanmean(pethMn,3);
subplot(nRow,2,1), hold on
ph=plot(pethBin,pethMnMn(:,nwIx));
set(ph,'linewidth',1.2)
niceyuax;
set(gca,'xlim',pethBin([1 end]));
xlabel('time (ms)');
ylabel('spx/bin');
% cum spx
subplot(nRow,2,2), hold on
ph=plot(pethBin,cumsum(pethMnMn(:,nwIx)));
set(ph,'linewidth',1.2)
niceyuax;
set(gca,'xlim',pethBin([1 end]));
xlabel('time (ms)');
ylabel('cumsum(spx/bin)');
legend(ph,num2str(indepParLevel_nw),'Location','Southeast');

% collection of handles for scalar results: experiments in rows, dependent
% parameters in columns
phArr=nan(nExpChanName,nDepPar);
% loop over scalar results
for k=1:nDepPar
  % normalized version? 
  if ap.depPar{k,2}  
    y=eval([ap.depPar{k,1} '_norm;']);
  else
    y=eval([ap.depPar{k,1} ';']);
  end
  % as we don't want any wash values in the summary plot delete them right
  % away
  y=y(:,nwIx,:);
  y=permute(y,[3 2 1]);
  subplot(nRow,nCol,k+scalarSpOffs), hold on
  finIx=isfinite(y);
  % set all infs to nan to avoid derailing of averages below
  y(~finIx)=nan;
  % ** set callback: culture ID to be displayed if mouse click on data point occurs
  % (mind the nans & infs - will give trouble when setting callbacks!)
  for g=1:nExpChanName
    ph=plot(indepParLevel_nw(finIx(g,:)),y(g,finIx(g,:)),'o-');
    if ~isempty(ph)
      phArr(g,k)=ph;
      set(ph,'color',[.8 .8 .8]);
      set(ph,'ButtonDownFcn',callB{g});
    end
  end
  % mean, std, 95% CI
  mn=nanmean(y,1);
  sd=nanstd(y,0,1);
  n=sum(finIx);
  ci=sd./sqrt(n).* -tinv(.05/2,n-1);
  if strcmp(ebType,'sd')
    ebh=errorbar(indepParLevel_nw',mn,sd,'ko-');
  else
    ebh=errorbar(indepParLevel_nw',mn,ci,'ko-');
  end    
  set(ebh','linewidth',1.6,'markersize',6,'markerfacecolor','k');
  hold off
  nicexyax(10);
  set(gca,'xtick',indepParLevel_nw);
  th=title(ap.depPar{k,1});
end

if ~isempty(ap.printFig)
  print(ap.printFig,'-r400',[ap.resPath ap.resFn '_SPX_sum']); 
end

% after printout set color
for g=1:nExpChanName
  set(phArr(g,isfinite(phArr(g,:))),'color',cm(g,:));
end
% ..and save as fig
saveas(gcf,[ap.resPath ap.resFn '_SPX_sum'],'fig');


% *************************************************************************
% *                 individual experiments' figures                       *
% *************************************************************************

if ap.IndividExperimentPlot
  % layout of figure:
  % - treat raw peths independent of scalar results
  % - reserve 
  % -- first row of figure for instfr 
  % -- second row for peth and cumpeth 
  if nDepPar==1
    nCol=1;
    nRow=3;
  elseif nDepPar<=6
    nCol=2;
    nRow=ceil(nDepPar/nCol)+2;
  else
    nCol=3;
    nRow=ceil(nDepPar/nCol)+2;
  end
  % the subplot offset for peth plots
  pethSpOffs=2;
  % the subplot offset for scalar plots
  scalarSpOffs=nCol*2;
  
  fh=figure(2); orient tall
  clf
  ftag='individual plot';
  labelscale('fontSz',8,'scaleFac',1.0,'lineW',1.2,'markSz',6);
  tmpScrSz=get(0,'Screensize');
  tmpScrSz([1 2])=tmpScrSz([1 2])+tmpScrSz([3 4]).*[.45 .02];
  tmpScrSz([3 4])=tmpScrSz([3 4]).*[.5 .9];
  set(fh,'position',tmpScrSz,'tag',ftag,'name',ftag,'color',[.9 .9 1],'numbertitle','off');
  
  for g=1:nExpChanName
%     % i) raster plots
%     for j=1:indepParNLevel
%       subplot(nRow,indepParNLevel,j), cla
%       % all in sec
%       tsl=Tsl{1,j,g}/1000;
%       % x interval of 20 s
%       x=mod(tsl,20);
%       % this operation implies a horizontal offset of 1 between consecutive rows
%       y= -floor(tsl/20);
%       ph=plot(x,y,'k.');
%       set(ph,'markersize',5);
%       niceyax(15);
%       set(gca,'xlim',[0 20],'ytick',[]);
%     end

    % i) instfr 
    subplot(nRow,1,1)
    cla, hold on
    % get standard colororder 
    colorOrd=get(gca,'colorOrder');
    ifrD={[]};
    dCt=0;
    for j=[nwIx washIx]
      if ~isempty(Tsl{1,j,g})
        dCt=dCt+1;
        % unit: s; 
        tsl=Tsl{1,j,g}/1000;
        ifrBin=0:ifrBinW:tsl(end);
        ifrD{dCt}=histc(tsl,ifrBin)/ifrBinW;
      end
    end
    % ** proceed only if dct~=0 
    if dCt
      pllplot(ifrD,'si',ifrBinW*1e6);
      %     ph=stairs(ifrBin,n);
      %     set(ph,'color',colorOrd(mod(j-1,size(colorOrd,1))+1,:));
      %     xlabel('time (s)');
      %     ylabel('Freq (Hz)')
      %     box on
      %     niceyuax;

      % ii) peth & cum peth
      subplot(nRow,2,pethSpOffs+1)
      ph=plot(pethBin,pethMn(:,[nwIx washIx],g));
      set(ph,'linewidth',1.2)
      niceyuax;
      set(gca,'xlim',pethBin([1 end]));
      xlabel('time (ms)');
      ylabel('spx/bin');
      % cum peth
      subplot(nRow,2,pethSpOffs+2)
      ph=plot(pethBin,cumsum(pethMn(:,[nwIx washIx],g)));
      set(ph,'linewidth',1.2)
      niceyuax;
      set(gca,'xlim',pethBin([1 end]));
      xlabel('time (ms)');
      ylabel('cumsum(spx/bin)');
      legend(ph,num2str(indepParLevel([nwIx washIx])),'Location','Southeast');

      % iii) loop over scalar results
      for k=1:nDepPar
        % plotting normalized data would not make a lot of sense
        y=eval([ap.depPar{k,1} ';']);
        y=y(1,:,g)';
        subplot(nRow,nCol,k+scalarSpOffs)
        % all except wash
        plot(indepParLevel(nwIx),y(nwIx),'ko-');
        % overplot point representing control value
        hold on
        ph=plot(indepParLevel(indepParNormIx),y(indepParNormIx),'ko-');
        set(ph,'MarkerFaceColor','k');
        if doWash
          % plot wash value at norm abscissa value as triangle
          plot(indepParLevel(indepParNormIx),y(washIx),'kv');
        end
        hold off
        nicexyax(8);
        % set lower y limit to zero so that changes are somwhat put into
        % perspective
        yl=get(gca,'ylim');
        set(gca,'ylim',[0 yl(2)]);
        set(gca,'xtick',indepParLevel);
        th=title(ap.depPar{k,1});
      end
      figtitle(expChanName{g});
      % print & 'refresh' figure, otherwise too much load
      if ~isempty(ap.printFig)
        % devise a friendly file name from expChanName: replace comma by
        % underscore and remove whitespace
        pfn=strrep(expChanName{g}(~isspace(expChanName{g})),',','_');
        pfn=[pfn '_SPX'];
        print(ap.printFig,'-r500',[ap.resPath pfn]);
      end
      drawnow
    else
      warning('skipping plot, no spx found for current experiment')
    end
  end
end

