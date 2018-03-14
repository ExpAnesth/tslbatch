function tslbatch_plotdata_fp(fn,varargin)
% ** function tslbatch_plotdata_fp(fn,varargin)
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
% (re-) compute some vars
% *** all FP scalar parameters
[tmp,dpIx]=union(ap.depPar(:,1),masterDepPar(cat(1,masterDepPar{:,3})==2,1));
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
% - target sampling freq of burst cutout figure (Hz) in individual plots
fs_buCutout=200;
% - length of burst cutout figure (s) in individual plots
len_buCutout=1;
% - create figure with return maps?
doReturnMap=false;

% *************************************************************************
% *                    summary figure                                     *
% *************************************************************************
indepParLevel_nw=indepParLevel(nwIx);
% --- layout of figure:
% - treat spectra independent of scalar results
% - reserve first row of figure for spectra
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

fh=figure(1); orient tall, clf
ftag='summary plot';
labelscale('fontSz',8,'scaleFac',1.0,'lineW',1.2,'markSz',6);
tmpScrSz=get(0,'Screensize');
tmpScrSz([1 2])=tmpScrSz([1 2])+tmpScrSz([3 4])*.02;
tmpScrSz([3 4])=tmpScrSz([3 4]).*[.5 .9];
set(fh,'position',tmpScrSz,'tag',ftag,'name',ftag,'color',[.9 .9 1],'numbertitle','off');
cm=coma('jet','ncols',nExpChanName);

% define simple callbacks for mouse click on idividual data points
callB=cell(nExpChanName,1);
for g=1:nExpChanName
  % callB{g}=['title(''' expChanName{g} ''',''interpreter'',''none'');'];
  callB{g}=['disp(''' expChanName{g} ''');'];
end

% plot psDens
psDensMnMn=nanmean(psDensMn,3);
subplot(nRow,2,1), hold on
ph=plot(psdF,psDensMnMn(:,nwIx));
set(ph,'linewidth',1.2)
set(gca,'xlim',psdF([1 end]),'xscale','log','yscale','lin');
% omit x axis label - will shrink axis
% xlabel('Freq (Hz)');
ylabel('psd');
legend(ph,num2str(indepParLevel_nw),'Location','eastoutside');

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

if ~isempty(ap.printFig), 
  print(ap.printFig,'-r400',[ap.resPath ap.resFn '_FP_sum']); 
end

% after printout set color
for g=1:nExpChanName
  set(phArr(g,isfinite(phArr(g,:))),'color',cm(g,:));
end
% ..and save as fig
saveas(gcf,[ap.resPath ap.resFn '_FP_sum'],'fig');


% *************************************************************************
% *                 individual experiments' figures                       *
% *************************************************************************

if ap.IndividExperimentPlot
  % layout of standard figure:
  % - treat spectra independent of scalar results
  % - reserve:
  % -- first row of figure for for raw plots 
  % -- second row for burst cutouts
  % -- third row for spectra
  
  % - semi-explicit assignment of rows and columns
  if nDepPar==1
    nCol=1;
    nRow=4;
  elseif nDepPar<=6
    nCol=2;
    nRow=ceil(nDepPar/nCol)+3;
  else
    nCol=3;
    nRow=ceil(nDepPar/nCol)+3;
  end
  % the subplot offset for burst cutout plot
  buExcSpOffs=2;
  % the subplot offset for psDens & buLen tc plot
  psDensSpOffs=4;
  % the subplot offset for scalar plots
  scalarSpOffs=nCol*3;
  
  fh_indiv=figure(2); orient tall, clf
  ftag='individual plot';
  labelscale('fontSz',8,'scaleFac',1.0,'lineW',1.2,'markSz',6);
  tmpScrSz=get(0,'Screensize');
  tmpScrSz([1 2])=tmpScrSz([1 2])+tmpScrSz([3 4]).*[.45 .02];
  tmpScrSz([3 4])=tmpScrSz([3 4]).*[.5 .9];
  set(fh_indiv,'position',tmpScrSz,'tag',ftag,'name',ftag,'color',[.9 .9 1],'numbertitle','off');
  
  if doReturnMap
    % figure for return maps
    fh_indivExt=figure(3); orient tall, clf
    ftag='individual plot - extended';
    labelscale('fontSz',8,'scaleFac',1.0,'lineW',1.2,'markSz',6);
    tmpScrSz=get(0,'Screensize');
    tmpScrSz([1 2])=tmpScrSz([1 2])+tmpScrSz([3 4]).*[.35 .02];
    tmpScrSz([3 4])=tmpScrSz([3 4]).*[.5 .9];
    set(fh_indivExt,'position',tmpScrSz,'tag',ftag,'name',ftag,'color',[.9 .9 1],'numbertitle','off');
  end
  
  for g=1:nExpChanName
    figure(fh_indiv);
    clf
    isRawDataFound=false;
    % container for complete recordings of current experiment
    dd={[]};
    % container for arrays of nans (as long as recording) with bursts
    % embedded
    dd_bu=dd;
    % container for bursts
    d_bu=dd;
    % index 
    ddIx=0;
    % wash recording last
    for j=[nwIx washIx]
      curRawFn=[fileInfo(1,j,g).dDir '\' fileInfo(1,j,g).fileName '.mat'];
      if exist(curRawFn,'file')
        isRawDataFound=true;
        ddIx=ddIx+1;
        tmpCh={expChanName{g}(strfind(expChanName{g},',')+2:end)};
        [d,si]=matDload(curRawFn,'channels',tmpCh);
        % downsample to ~fs_buCutout Hz for these coarse overview plots, no matter what
        sampFac=floor((1e6/fs_buCutout)/si);
        si=si*sampFac;
        d=d(1:sampFac:end);
        dd{ddIx}=d;
        dd_bu{ddIx}=repmat(nan,size(d));
        % index to cutouts:
        % - exact
        tmpBuIx=cumsum(Etsl{1,j,g}(:,[etslc.tsCol etslc.durCol]),2);
        buIx=cont2discrete(tmpBuIx,si/1000);
        % - shifted by 50 ms to left so that base line can be seen
        buIx_shift=cont2discrete(tmpBuIx-50,si/1000);        
        % preallocate (1 s length)
        len_buCutout_pts=floor(len_buCutout/(si/1e6));
        d_buExc{ddIx}=repmat(nan,len_buCutout_pts,size(buIx,1));
        for cix=1:size(buIx,1)
          % embed
          dd_bu{ddIx}(buIx(cix,1):buIx(cix,2))=d(buIx(cix,1):buIx(cix,2));
          % collect up to length of 1 s (shifted version)
          redIx=1:min(diff(buIx_shift(cix,[1 2])),len_buCutout_pts);
          d_buExc{ddIx}(redIx,cix)=d(buIx_shift(cix,1):buIx_shift(cix,1)+redIx(end)-1);
        end
      end
    end

    % i) plots of bursts embedded in recording
    sph=subplot(nRow,1,1);
    cla
    if isRawDataFound
      [ylim,dy,yscaleFac,ph]=pllplot(dd,'si',si,'noscb',1);
      set(ph,'linewidth',.25);
      hold on
      [ylim,dy,yscaleFac,ph]=pllplot(dd_bu,'si',si,'ylim',ylim,'spacing','fixed','dy',dy,'noscb',1);
      set(ph,'linewidth',.3,'color','m');
    end

    % ii) burst excerpts
    if isRawDataFound
      axhArr=repmat(nan,1,ddIx);
      for j=1:ddIx
        axhArr(j)=subplot(nRow,ddIx,j+ddIx);
        if ~isempty(d_buExc{j})
          [ylim,dy,yscaleFac,ph]=pllplot(d_buExc{j},'si',si,'spacing','fixed','dy',.000001,'noscb',1);
          set(ph,'linewidth',.35);
          axis tight
        end
      end
      subpax(gcf,'spInd',axhArr);
    end
    
    % iii) time course of bu len
    subplot(nRow,2,1+psDensSpOffs)
    % get standard colororder 
    colorOrd=get(gca,'colorOrder');
    cla, hold on
    ct=0;
    for j=[nwIx washIx]
      ct=ct+1;
      etsl=Etsl{1,j,g};
      ph=plot(etsl(:,etslc.tsCol)/1000,etsl(:,etslc.durCol)/1000);
      set(ph,'linewidth',1.2,'color',colorOrd(ct,:),'marker','v','markersize',4);
    end
    niceyuax;
    xlabel('time (s)');
    ylabel('burst len (s)');

    % iv) psDens
    subplot(nRow,2,2+psDensSpOffs)
    cla, hold on
    % get standard colororder 
    colorOrd=get(gca,'colorOrder');
    ph=plot(psdF,psDensMn(:,[nwIx washIx],g));
    set(ph,'linewidth',1.2)
    niceyuax;
    set(gca,'xlim',psdF([1 end]),'xscale','log','yscale','lin');
    xlabel('Freq (Hz)');
    ylabel('psd');
    legend(ph,num2str(indepParLevel([nwIx washIx])),'Location','EastOutside');

    % v) loop over scalar results
    for k=1:nDepPar
      % plotting normalized data would not make a lot of sense
      y=eval([ap.depPar{k,1} ';']);
      y=y(1,:,g)';
      subplot(nRow,nCol,k+scalarSpOffs)
      % all except wash
      ph=plot(indepParLevel(nwIx),y(nwIx),'ko-');
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
    if ~isempty(ap.printFig),
      % devise a friendly file name from expChanName: replace comma by
      % underscore and remove whitespace
      pfn=strrep(expChanName{g}(~isspace(expChanName{g})),',','_');
      pfn=[pfn '_FP'];
      print(ap.printFig,'-r500',[ap.resPath pfn]);
    end
    drawnow
    % ---------- return maps (optional) -----------------------------------
    if doReturnMap
      figure(fh_indivExt); 
      clf
      % - length silent period vs. length preceding burst 
      subplot(2,2,1)
      hold on
      for h=1:indepParNLevel
        etsl=Etsl{1,h,g};
        nBu=size(etsl,1);
        silentEtsl=SilentEtsl{1,h,g};
        nSp=size(silentEtsl,1);
        plotArr=repmat(nan,max(nBu,nSp),2);
        pCol=colorOrd(mod(h-1,size(colorOrd,1))+1,:);
        if nSp>1 && nBu>0
          plotArr(1:nBu,1)=etsl(:,etslc.durCol);
          % if first burst comes before first silent period...
          if etsl(1,etslc.tsCol)<silentEtsl(1,etslc.tsCol)
            plotArr(1:nSp,2)=silentEtsl(:,etslc.durCol);
          else
            plotArr(1:nSp-1,2)=silentEtsl(2:end,etslc.durCol);
          end
          ph=plot(plotArr(:,1),plotArr(:,2),'o');
          set(ph,'color',pCol,'markerfacecolor',pCol,'markersize',4);
        end
      end
      nicexyax;
      xlabel('length bu(n-1)');
      ylabel('length sp(n)');
      hold off

      % - length burst vs. length preceding burst
      subplot(2,2,2)
      hold on
      for h=1:indepParNLevel
        etsl=Etsl{1,h,g};
        nBu=size(etsl,1);
        pCol=colorOrd(mod(h-1,size(colorOrd,1))+1,:);        
        if nBu>1
          ph=plot(etsl(2:end,etslc.durCol),etsl(1:end-1,etslc.durCol),'o');
          set(ph,'color',pCol,'markerfacecolor',pCol,'markersize',4);
        end
      end
      nicexyax;
      xlabel('length bu(n-1)');
      ylabel('length bu(n)');
      hold off
      % §§ legend, print
    end    
  end
end


