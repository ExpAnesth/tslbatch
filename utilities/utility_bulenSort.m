clear


% mat file name (without extension)
mFn='d:\hh\_tbs\AChpilot.mat';
mFn='d:\hh\projects\ctx_gabaUptake\resultsnfigs_SNAP\tslbatchPlots\tslbatch_SNAP.mat';

% levels of independent variable (=concentrations in most cases) to deal
% with IN ASCENDING ORDER
curIpVal=[0];

printas=[];  % don't save graphics
% printas='-djpeg97'; % save as jpeg ('working style')
% printas='-dpsc2'; % save as color post script (nice finish)

% graphics options
labelscale('fontSz',10,'scaleFac',1,'lineW',.25,'markSz',4);
ornt='tall';
% name & full path to figure to be saved
figName=['D:\nix'];

etslconst;
nIpVal=numel(curIpVal);
% load *.mat data file
load(mFn);

% in most cases, the independent parameter is the concentration of the
% substance in question. This command plots its levels of on screen
indepParLevel
% index to columns in all variables
[nix,columnIx]=intersect(indepParLevel,curIpVal);
if ~isequal(nix,curIpVal)
  errordlg({'variable ''curIpVal'' contains values that do not exist in the data.',...
    'Available values: ', num2str(indepParLevel)})
  return
end

Etsl=permute(Etsl(:,columnIx,:),[3 2 1]);
mdBurstLen=permute(mdBurstLen(:,columnIx,:),[3 2 1])/1000;
% maxBurstLen=permute(maxBurstLen(:,columnIx,:),[3 2 1])/1000;
nExp=size(Etsl,1);

figure(2), clf;
plotType='scatter';

switch plotType
  case 'scatter'
    % - each burst by a circle
    % - median 
    % - interquartile range by thick solid line
    % - plot on abscissa according to median 
    % - histogram
    sph=subplot('position',[.35 .05 .6 .9]);
    hold on
    catBuLen=[];
    for g=1:nExp
      buLen=Etsl{g}(:,etslc.durCol)/1000;
      catBuLen=cat(1,catBuLen,buLen);
      if ~isempty(buLen)
        iqr=prctile(buLen,[25 75]);
        % §§
        maBuLe=max(buLen);
        lh=line(mdBurstLen(g)*[1 1],iqr,'color','k');
        set(lh,'linewidth',2)
        plot(mdBurstLen(g),buLen,'ko');
        % max burst len
        ph=plot(mdBurstLen(g),maBuLe,'bo');
        set(ph,'markerfacecolor','b');
        % median burst len
        ph=plot(mdBurstLen(g),mdBurstLen(g),'ko');
        set(ph,'markerfacecolor','r','markersize',7);
      end
    end
    set(gca,'xscale','log');
    nicexyax(40);
    yl=get(gca,'ylim');
    xlabel('median burst length (s)')
    ylabel('burst length (s)')
    % histogram
    sph=subplot('position',[.05 .05 .25 .9]);
    hold on
    binW=.5;
    bin=0:binW:yl(2);
    % distrib of median bu len
    nMd=histc(mdBurstLen,bin);
    barh(bin+binW/2,nMd,1.0,'k')
%     % distrib of bu len (all recordings)
%     n=histc(catBuLen,bin);
%     contourbarh(bin,n,'color','r');
    set(gca,'ylim',yl,'xdir','reverse');
    
    
  case 'distrPlot'
    collBuLen=cell(1,nExp);
    for g=1:nExp
      buLen=Etsl{g}(:,etslc.durCol)/1000;
      collBuLen{g}=buLen;
    end
    % sort according to median
    [mdBurstLen,ix]=sort(mdBurstLen);
    collBuLen=collBuLen(ix);
    % plain hist with fixed bins
    hand=distributionPlot(collBuLen,1.7,3,[],0,[.25:.5:50]);
    %   % smoothed
    %   hand=distributionPlot(collBuLen,1.9,3,[],1,10);
    % set(hand{,
    
  case 'boxPlot'
    catBuLen=[];
    group=[];
    catMdBuLen=[];
    for g=1:nExp
      buLen=Etsl{g}(:,etslc.durCol);
      catBuLen=cat(1,catBuLen,buLen);
      group=cat(1,group,repmat(g,size(buLen,1),1));
      catMdBuLen=cat(1,catMdBuLen,repmat(mdBurstLen(g),size(buLen,1),1));
    end
    % sort according to median
    [nix,ix]=sort(catMdBuLen);
    
    boxplot(catBuLen,-catMdBuLen);
    % set(gca,'xtick',1:nExp,'xticklabel',int2str(flipud(unique(group))))
    
    % % as boxplot automatically sorts values according to values in 'group', ...
    % boxplot(catBuLen(ix,:),group(ix))
    %
    
end

% ----- printing -----
if ~isempty(printas), 
  print(printas,figName,'-r400'); 
end