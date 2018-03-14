function tslbatch_timecourse_fp(fn,varargin)
% ** function tslbatch_timecourse_fp(fn,varargin)

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
nDepPar=size(ap.depPar,1);
nExpChanName=numel(expChanName);
indepParNLevel=numel(indepParLevel);

etslconst;

labelscale('fontSz',10,'scaleFac',1.0,'lineW',.5,'markSz',5);

for g=1:nExpChanName
  etsl_collect=[];
  silentEtsl_collect=[];
  % time in seconds from midnight at which culture was exposed to ACSF -
  % this information should be in undef2. No safeguard whatsoever against
  % nonsense values
  tmp=floor(undef2(1,1,g)/100);
  tOffset=tmp*60*60 + 60*(undef2(1,1,g)-tmp*100);
  figure(g), clf, orient landscape
  axh1=gca;
  axh2=axes('pos',get(axh1,'pos'));
  set([axh1 axh2],'nextplot','add');
  set(axh2,'yaxisloc','right','xtick',[],'color','none');
  for fIx=1:indepParNLevel
    etsl=Etsl{1,fIx,g};
    silentEtsl=SilentEtsl{1,fIx,g};
    % string for callback
    cb=['disp(''' [fileInfo(1,fIx,g).fileName ', ' fileInfo(1,fIx,g).chName] ''');'];
    if ~isempty(etsl)
      etsl(:,etslc.tsCol)=etsl(:,etslc.tsCol)+(fileInfo(1,fIx,g).recTime(1)-tOffset)*1000;
      % insert row of nans to have disjunct recordings disconnected in plot
      etsl(end+1,:)=nan;
      axes(axh1);
      h1=plot(etsl(:,etslc.tsCol)/1000,etsl(:,etslc.durCol),'b-o');
      set(h1,'ButtonDownFcn',cb);
    end
    if ~isempty(silentEtsl)
      silentEtsl(:,etslc.tsCol)=silentEtsl(:,etslc.tsCol)+(fileInfo(1,fIx,g).recTime(1)-tOffset)*1000;
      silentEtsl(end+1,:)=nan;
      axes(axh2);
      h2=plot(silentEtsl(:,etslc.tsCol)/1000,silentEtsl(:,etslc.durCol),'g-^');
      set(h2,'ButtonDownFcn',cb);
    end
  end
  
  title(expChanName{g},'interpreter','none');
  
  axes(axh1);
  axis tight
  xlabel('time (s)')
  ylabel('burst length (ms)');
  axes(axh2);
  axis tight
  ylabel('silent period length (ms)');

  %xl1=get(axh1,'xlim');
  %xl2=get(axh2,'xlim');
  %xl=[0 max([xl1 xl2])];
  %set([axh1 axh2],'xlim',xl);

  if ~isempty(ap.printFig),
    print(ap.printFig,[ap.resPath 'tc_' expChanName{g}]);
  end
end



