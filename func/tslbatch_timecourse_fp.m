function tslbatch_timecourse_fp(fn,varargin)
% ** function tslbatch_timecourse_fp(fn,varargin)
% is a relatively quick-and-dirty function intended to plot the development
% over time (across many files) of some key parameters of neuronal
% activity. It plots LFP and SPX based parameters depending on
% availability. 
% IMPORTANT: in the file setting up tslbatch parameters, ap.indepPar must
% NOT be drug concentration but drug application order as in this examples:
%
% ap.indepPar={...
%   'drug1_applicRank',0
%   };
% 
% Furhermore, before calling tslbatch_timecourse_fp you must call
% tslbatch_procdata_fp 


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

labelscale('fontSz',10,'scaleFac',1.0,'lineW',.5,'markSz',3);

% 'bin width' (=interval for calculation) for firing rates and ISI
% computation in s
binW=60;
% bins of logarithms of ISI, capturing all values between 148 and 22000 ms
binLogISI=5:.1:10;


for g=1:nExpChanName
  % time in seconds from midnight at which culture was exposed to ACSF -
  % this information should be in undef2. No safeguard whatsoever against
  % nonsense values
  tmp=floor(undef2(1,1,g)/100);
  tOffset=tmp*60*60 + 60*(undef2(1,1,g)-tmp*100);
  figure(g), clf, orient landscape
  colormap(flipud(gray));
  
  % first subplot for firing rate and relative active time
  sph1=subplot(2,1,1);
  fixIntvAxH1=sph1;
  set(fixIntvAxH1,'color','none')
  fixIntvAxH2=axes('pos',get(fixIntvAxH1,'pos'));
  set([fixIntvAxH1 fixIntvAxH2],'nextplot','add');
  set(fixIntvAxH2,'yaxisloc','right','xtick',[],'color','none');

  % second subplot for burst and silent period durations
  sph2=subplot(2,1,2);
  buAxH1=sph2;
  buAxH2=axes('pos',get(buAxH1,'pos'));
  set([buAxH1 buAxH2],'nextplot','add');
  set(buAxH2,'yaxisloc','right','xtick',[],'color','none');
  
  
  % all (silent) etsls, concatenated
  etsl_cat=[];
  silentEtsl_cat=[]; 
  % spx tsl
  tsl_cat=[];
  for fIx=1:indepParNLevel
    etsl=Etsl{1,fIx,g};
    silentEtsl=SilentEtsl{1,fIx,g};
    tsl=Tsl{1,fIx,g};
    % string for callback
    cb=['disp(''' [fileInfo(1,fIx,g).fileName ', ' fileInfo(1,fIx,g).chName] ''');'];
    % --- SPX
    if ~isempty(tsl)
      % in keeping with code below, which was written earlier, compute &
      % plot on the spot here as well:
      bin=0:binW:round(diff(fileInfo(1,fIx,g).recTime));
      
      % - average firing freq in Hz
      hst=histc(tsl/1000,bin)/binW;
      
      % - ISI
      isiH=nan(numel(binLogISI),numel(bin)-1);
      for ii=1:numel(bin)-1
        partTsl=tsl(tsl>=bin(ii)*1000 & tsl<bin(ii+1)*1000);
        isiH(:,ii)=histc(log(diff(partTsl)),binLogISI);
      end
      
      % don't forget to kill last bin of firing rate hist
      hst(end)=[];
      bin(end)=[];
      
      % abscissa values
      absci=bin+binW/2+fileInfo(1,fIx,g).recTime(1)-tOffset;
      
      % plot into subplot with bin centers as abscissa
      axes(fixIntvAxH2);
      imagesc(absci,binLogISI,sqrt(isiH));
      axes(fixIntvAxH1);
      ph=plot(absci,hst,'k-o');
      set(ph,'markersize',5,'markerfacecolor','r');
      
      % plot value of indep par (conc), color-coded
      ph=plot(absci,0,'ks');
      set(ph,'markersize',8,'markerfacecolor',[0 1-drug1_conc(1,fIx,g)/100 0]);
           
      % ** automatically determine bursts and silent periods from spx if no
      % genuine burst data available
      if isempty(etsl) && isempty(silentEtsl)
        [etsl,~,silentEtsl,stats]=etslburstf(tsl,200,'maxIEI_init',50,...
          'minNEvPerBurst',2,'recLen',diff(fileInfo(1,fIx,g).recTime)*1000);
%         % relative time in burst
%         axes(fixIntvAxH2);
%         ph=plot(mean(bin)+fileInfo(1,fIx,g).recTime(1)-tOffset,stats.relTimeInBurst,'kd');
%         set(ph,'markersize',8,'markerfacecolor','c');
      end
      
    end
    
    if ~isempty(silentEtsl)
      silentEtsl(:,etslc.tsCol)=silentEtsl(:,etslc.tsCol)+(fileInfo(1,fIx,g).recTime(1)-tOffset)*1000;
      silentEtsl(end+1,:)=nan;
      silentEtsl_cat=cat(1,silentEtsl_cat,silentEtsl);
      axes(buAxH1);
      % § plot in gray so that burst length is more prominent
      h2=plot(silentEtsl(:,etslc.tsCol)/1000,conv(silentEtsl(:,etslc.durCol),triang(9),'same')/sum(triang(9)),'w-^');
      set(h2,'ButtonDownFcn',cb,'color',[.7 .7 .7]);
    end
    
    if ~isempty(etsl)
      etsl(:,etslc.tsCol)=etsl(:,etslc.tsCol)+(fileInfo(1,fIx,g).recTime(1)-tOffset)*1000;
      % insert row of nans to have disjunct recordings disconnected in plot
      etsl(end+1,:)=nan;
      etsl_cat=cat(1,etsl_cat,etsl);
      axes(buAxH2);
      % ** smoothed burst length!
      h1=plot(etsl(:,etslc.tsCol)/1000,conv(etsl(:,etslc.durCol),triang(9),'same')/sum(triang(9)),'b-o');
      set(h1,'ButtonDownFcn',cb);
    end
    
    
  end
  
  axes(sph1);
  title(expChanName{g},'interpreter','none');
  
  % plot parameters based on fixed interval analysis
  axes(fixIntvAxH1)  
  ylabel('average firing rate');
  nicexy0ax
  grid on
  % common x axis for all plots
  xl=get(fixIntvAxH1,'xlim');
  xl=[0 xl(2)];
  xti=0:600:6e4;
  set(fixIntvAxH1,'xlim',xl,'xtick',xti);
  axes(fixIntvAxH2)
  axis tight
  set(fixIntvAxH2,'xlim',xl,'xtick',[]);
  set(fixIntvAxH2,'ytick',log([250 1000 4000 16000]),'yticklabel',[.25 1 4 16]);
  ylabel('ISI (s)');
  % focus back on first axis
  axes(fixIntvAxH1)
  
%   nicexy0ax
%   set(fixIntvAxH2,'xlim',get(buAxH1,'xlim'),'xtick',get(buAxH1,'xtick'));
%   ylabel('rel act time');
  
  % compute percentiles for y axis of burst/sp duration so as to exclude
  % extremes from scaling
  yl_etsl=prctile(etsl_cat(:,etslc.durCol),99);
  yl_silentEtsl=prctile(silentEtsl_cat(:,etslc.durCol),99);  
  axes(buAxH2);
  axis tight
  set(buAxH2,'ylim',[0 yl_etsl],'xlim',xl,'xtick',xti);
  xlabel('time (s)')
  ylabel('burst length (ms)');
  grid on
  axes(buAxH1);
  axis tight
  set(buAxH1,'ylim',[0 yl_silentEtsl],'xlim',xl,'xtick',[]);
  ylabel('silent period length (ms)');
  % focus back on second (transparent) axis
  axes(buAxH2);
  
  
  
  if ~isempty(ap.printFig),
    print(ap.printFig,[ap.resPath 'tc_' expChanName{g}]);
  end
end



