function tslbatch_plotdata(fn,varargin)
% ** function tslbatch_plotdata(fn,varargin)
% plotting routine of tslbatch workhorse. See the documentation.

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
[tmp,dpIx]=intersect(ap.depPar(:,1),masterDepPar(cat(1,masterDepPar{:,3})==2,1));
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
% background color of axes and figures
axBgcolor=[.88 .88 .88];
figBgcolor=[.85 .85 .85];
% factors by which to blow up axes
xBlowupFac=1.05;
yBlowupFac=1.05;

% *************************************************************************
% *                    summary figure                                     *
% *************************************************************************
indepParLevel_nw=indepParLevel(nwIx);
% --- layout of figure:
% - treat raw peths independent of scalar results
% - reserve first row of figure for peths
% - semi-explicit assignment of rows and columns
if nDepPar<=6
  nCol=3;
  nRow=3;
elseif nDepPar<=8
  nCol=4;
  nRow=ceil(nDepPar/nCol)+1;
else
  nCol=5;
  nRow=ceil(nDepPar/nCol)+1;
end
% the subplot offset for scalar plots
scalarSpOffs=nCol;

fh=figure(1); orient landscape, clf
ftag='summary plot';
labelscale('fontSz',8,'scaleFac',1.0,'lineW',1.2,'markSz',6);
tmpScrSz=get(0,'Screensize');
% tmpScrSz([1 2])=tmpScrSz([1 2])+tmpScrSz([3 4])*.02;
% tmpScrSz([3 4])=tmpScrSz([3 4]).*[.5 .9];
tmpScrSz([1 2])=tmpScrSz([1 2])+tmpScrSz([3 4]).*[.02 .25];
tmpScrSz([3 4])=tmpScrSz([3 4]).*[.9 .65];
set(fh,'position',tmpScrSz,'tag',ftag,'name',ftag,'color',figBgcolor,'numbertitle','off');
% colormap for plots in which color codes for individual experiments
cm=coma('jet','ncols',nExpChanName);
% colormap for plots in which color codes for independent par (conc mostly)
cmIdep=coma('jet','ncols',numel(nwIx));
if doWash
  % append gray for wash
  cmIdep=cat(1,cmIdep,[.6 .6 .6]);
end
% define simple callbacks for mouse click on idividual data points
callB=cell(nExpChanName,1);
for g=1:nExpChanName
  % callB{g}=['title(''' expChanName{g} ''',''interpreter'',''none'');'];
  callB{g}=['disp(''' expChanName{g} ''');'];
end

% plot peths:
% first, normalize pethN so that it reflects the proportion of bursts which
% are 'on' at each bin
pethN=pethN./repmat(max(pethN),[size(pethN,1) 1]);
% i. burst start
pethMnMn=nanmean(pethMn,3);
subplot(nRow,4,1), cla, hold on
rexy('ax',gca,'xfac',xBlowupFac,'yfac',yBlowupFac);
set(gca,'colororder',cmIdep,'color',axBgcolor);
% binW=diff(pethBin(1:2));
% ph=plot(pethBin,pethMnMn(:,nwIx)/binW*1000);
ph=plot(pethBin,pethMnMn(:,nwIx));
set(ph,'linewidth',1.4)
niceyuax(30)
set(gca,'xlim',pethBin([1 end]));
xlabel('time (ms)');
ylabel('firing rate (Hz)');
% legend(ph,num2str(indepParLevel_nw),'Location','Northeast');
for ii=1:numel(indepParLevel_nw)
  smarttext(num2str(indepParLevel_nw(ii)),.8,1-.09*ii,'fontweight','bold',...
    'fontsize',8,'HorizontalAlignment','right','color',cmIdep(ii,:));
end
title('peri-burst start')

% ii. burst prevalence during PE interval
pethNMn=nanmean(pethN,3);
subplot(nRow,4,2), cla, hold on
rexy('ax',gca,'xfac',xBlowupFac,'yfac',yBlowupFac);
set(gca,'colororder',cmIdep,'color',axBgcolor);
ph=plot(pethBin,pethNMn(:,nwIx));
set(ph,'linewidth',1.4)
niceyuax(30)
set(gca,'xlim',pethBin([1 end]));
xlabel('time (ms)');
ylabel('burst prevalence');
title('peri-burst start')

% iii. cum spx
subplot(nRow,4,3), cla, hold on
rexy('ax',gca,'xfac',xBlowupFac,'yfac',yBlowupFac);
set(gca,'colororder',cmIdep,'color',axBgcolor);
% ph=plot(pethBin,cumsum(pethMnMn(:,nwIx))/binW*1000);
ph=plot(pethBin,cumsum(pethMnMn(:,nwIx)));
set(ph,'linewidth',1.4)
niceyuax(30)
% y ticks useless with cum hist
set(gca,'xlim',pethBin([1 end]),'ytick',[]);
xlabel('time (ms)');
ylabel('cumsum(firing rate)');
title('peri-burst start')

% iv. peths of silent periods (=burst ends)
spPethMnMn=nanmean(spPethMn,3);
subplot(nRow,4,4), cla, hold on
rexy('ax',gca,'xfac',xBlowupFac,'yfac',yBlowupFac);
set(gca,'colororder',cmIdep,'color',axBgcolor);
% binW=diff(spPethBin(1:2));
% ph=plot(spPethBin,spPethMnMn(:,nwIx)/binW*1000);
ph=plot(spPethBin,spPethMnMn(:,nwIx));
set(ph,'linewidth',1.4)
niceyuax(30)
set(gca,'xlim',spPethBin([1 end]));
xlabel('time (ms)');
ylabel('firing rate (Hz)');
title('peri-burst end')

% collection of handles for scalar results: experiments in rows, dependent
% parameters in columns
phArr=repmat(nan,nExpChanName,nDepPar);
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
  rexy('ax',gca,'xfac',xBlowupFac,'yfac',yBlowupFac);
  set(gca,'color',axBgcolor);
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
  nicexyax(15);
  set(gca,'xtick',indepParLevel_nw);
  th=title(ap.depPar{k,1});
end

if ~isempty(ap.printFig), 
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
  % - reserve: 
  % -- first row of figure for instfr 
  % -- second row for peth and cumpeth 
  % -- third row for effect size of peth §§ plus maybe other things to come
  if nDepPar==1
    nCol=2;
    nRow=4;
  elseif nDepPar<=3
    nCol=3;
    nRow=ceil(nDepPar/nCol)+3;
  else
    nCol=4;
    nRow=ceil(nDepPar/nCol)+3;
  end
  % the subplot offset for peth plots
  pethSpOffs=4;
  % the subplot offset for effect size plots
  esSpOffs=6;
  % the subplot offset for scalar plots
  scalarSpOffs=nCol*3;
  
  fh=figure(2); orient landscape, clf
  ftag='individual plot';
  labelscale('fontSz',8,'scaleFac',1.0,'lineW',1.2,'markSz',6);
  tmpScrSz=get(0,'Screensize');
  % tmpScrSz([1 2])=tmpScrSz([1 2])+tmpScrSz([3 4]).*[.45 .02];
  % tmpScrSz([3 4])=tmpScrSz([3 4]).*[.5 .9];
  tmpScrSz([1 2])=tmpScrSz([1 2])+tmpScrSz([3 4]).*[.02 .1];
  tmpScrSz([3 4])=tmpScrSz([3 4]).*[.9 .65];
  set(fh,'position',tmpScrSz,'tag',ftag,'name',ftag,'color',figBgcolor,'numbertitle','off');
  
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

    % row 1: instfr 
    subplot(nRow,1,1)
    cla, hold on
    rexy('ax',gca,'xfac',1.1,'yfac',yBlowupFac);
    ifrD={[]};
    ct=0;
    for j=[nwIx washIx]
      if ~isempty(Tsl{1,j,g})
        ct=ct+1;
        % unit: s; 
        tsl=Tsl{1,j,g}/1000;
        ifrBin=0:ifrBinW:tsl(end);
        ifrD{ct}=histc(tsl,ifrBin)/ifrBinW;
      end
    end
    % ** proceed only if dct~=0 
    if ct
      % row1: inst fr
      [nada1,nada2,nada3,ph]=pllplot(ifrD,'si',ifrBinW*1e6);
      set(ph,'linewidth',.8);

      % row2: i) peth 
      subplot(nRow,4,pethSpOffs+1)
      cla, hold on
      rexy('ax',gca,'xfac',xBlowupFac,'yfac',yBlowupFac);
      set(gca,'colororder',cmIdep,'color',axBgcolor);
%       binW=diff(pethBin(1:2));
%       ph=plot(pethBin,pethMn(:,[nwIx washIx],g)/binW*1000);
      ph=plot(pethBin,pethMn(:,[nwIx washIx],g));
      set(ph,'linewidth',1.4)
      niceyuax(30)
      set(gca,'xlim',pethBin([1 end]));
      xlabel('time (ms)');
      ylabel('firing rate (Hz)');
      tmp=indepParLevel([nwIx washIx]);
      for ii=1:numel(tmp)
        smarttext(num2str(tmp(ii)),.8,1-.09*ii,'fontweight','bold',...
          'fontsize',8,'HorizontalAlignment','right','color',cmIdep(ii,:));
      end
      % legend(ph,num2str(indepParLevel([nwIx washIx])),'Location','Southeast');
      
      % row2: ii) burst prevalence
      subplot(nRow,4,pethSpOffs+2)
      cla, hold on
      rexy('ax',gca,'xfac',xBlowupFac,'yfac',yBlowupFac);
      set(gca,'colororder',cmIdep,'color',axBgcolor);
      ph=plot(pethBin,pethN(:,[nwIx washIx],g));
      set(ph,'linewidth',1.4)
      niceyuax(30)
      set(gca,'xlim',pethBin([1 end]));
      xlabel('time (ms)');
      ylabel('burst prevalence');
            
      % row2: iii) cum peth 
      subplot(nRow,4,pethSpOffs+3)
      cla, hold on
      rexy('ax',gca,'xfac',xBlowupFac,'yfac',yBlowupFac);
      set(gca,'colororder',cmIdep,'color',axBgcolor);
      % ph=plot(pethBin,cumsum(pethMn(:,[nwIx washIx],g))/binW*1000);
      ph=plot(pethBin,cumsum(pethMn(:,[nwIx washIx],g)));
      set(ph,'linewidth',1.4)
      niceyuax(30)
      set(gca,'xlim',pethBin([1 end]),'ytick',[]);
      xlabel('time (ms)');
      ylabel('cumsum(fr)');
      
      % row2: iv) peths of silent periods
      subplot(nRow,4,pethSpOffs+4)
      cla, hold on
      rexy('ax',gca,'xfac',xBlowupFac,'yfac',yBlowupFac);
      set(gca,'colororder',cmIdep,'color',axBgcolor);
      % binW=diff(spPethBin(1:2));
      % ph=plot(spPethBin,spPethMn(:,[nwIx washIx],g)/binW*1000);
      ph=plot(spPethBin,spPethMn(:,[nwIx washIx],g));
      set(ph,'linewidth',1.4)
      niceyuax(30)
      set(gca,'xlim',spPethBin([1 end]));
      xlabel('time (ms)');
      ylabel('firing rate (Hz)');

      % row3: i) effect size
      subplot(nRow,3,esSpOffs+1)
      cla, hold on
      rexy('ax',gca,'xfac',xBlowupFac,'yfac',yBlowupFac);
      set(gca,'colororder',cmIdep,'color',axBgcolor);
      esColl=repmat(nan,size(pethBin,1),numel([nwIx washIx]));
      ct=0;
      for j=[nwIx washIx]
        ct=ct+1;
        if j~=indepParNormIx && size(RawPeth{1,j,g},2)>1 && size(RawPeth{1,indepParNormIx,g},2)>1
          stats=mes(RawPeth{1,j,g}',RawPeth{1,indepParNormIx,g}','U3','nBoot',0);
          esColl(:,ct)=stats.U3';
        end
      end
      ph=plot(pethBin,esColl,'o-');
      set(ph,'linewidth',.5,'markersize',3)
      axis tight
      set(gca,'xlim',pethBin([1 end]));
      set(gca,'ygrid','on');
      lh=line(pethBin([1 end]),[0 0],'color','k','linewidth',1.4,'linestyle','--');
      xlabel('time (ms)');
      ylabel('ES');

      % row3: ii) CV of bin counts
      subplot(nRow,3,esSpOffs+2)
      cla, hold on
      rexy('ax',gca,'xfac',xBlowupFac,'yfac',yBlowupFac);
      set(gca,'colororder',cmIdep,'color',axBgcolor);
      % re-use variables from above
      esColl=repmat(nan,size(pethBin,1),numel([nwIx washIx]));
      ct=0;
      for j=[nwIx washIx]
        ct=ct+1;
        if ~isempty(RawPeth{1,j,g})
          esColl(:,ct)=nanstd(RawPeth{1,j,g},1,2)./nanmean(RawPeth{1,j,g},2);
        end
      end
      ph=plot(pethBin,esColl,'v-');
      set(ph,'linewidth',.5,'markersize',3)
      axis tight
      set(gca,'xlim',pethBin([1 end]));
      set(gca,'xgrid','on');
      xlabel('time (ms)');
      ylabel('CV bin count');

      % rows 4 and following: loop over scalar results
      for k=1:nDepPar
        % plotting normalized data would not make a lot of sense
        y=eval([ap.depPar{k,1} ';']);
        y=y(1,:,g)';
        subplot(nRow,nCol,k+scalarSpOffs)
        rexy('ax',gca,'xfac',xBlowupFac,'yfac',yBlowupFac);
        % all except wash
        ph=plot(indepParLevel(nwIx),y(nwIx),'ko-');
        % overplot point representing control value
        hold on
        ph=plot(indepParLevel(indepParNormIx),y(indepParNormIx),'ko-');
        set(ph,'MarkerFaceColor','k');
        set(gca,'color',axBgcolor);
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
      if ~isempty(ap.printFig),
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

