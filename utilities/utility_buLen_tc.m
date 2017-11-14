% plots time course of burst length for a given set of data
dStr='2009_07_14';
ddir=['\\Wks1\ACSF 3.0\' dStr '_SET2\'];
tOffset=9*60*60 + 27*60;
fIndex=[0 1 6 7 12 13 18];

dStr='2009_07_14';
ddir=['\\Wks1\ACSF 3.0\' dStr '_SET2\'];
fIndex=[20 21 26 27 32 33 38];
tOffset=11*60*60;

dStr='2009_07_17';
ddir=['\\Wks1\ACSF 3.0\' dStr '_SET2\'];
fIndex=[0 1 6 7 12 13 19];
tOffset=11*60*60;

dStr='2009_07_17';
ddir=['\\Wks1\ACSF 3.0\' dStr '_SET2\'];
fIndex=[22 23 28 29 34 35 40];
tOffset=12*60*60 + 25*60;

% by default, the beginning of the first recording is time zero; if you
% want it to have a different value, add it here (unit: s) 

numOfFiles=length(fIndex);
fn=[];
% from the vector of integers make a vector of strings with width 4, namely
% the file indices (zero-padding)
for i=1:numOfFiles
  fn(i,:)=[dStr '_' sprintf('%.4i',fIndex(i))  '_IN7_FP_res.mat']; 
end

% graphics options
labelscale('fontSz',10,'scaleFac',1,'lineW',.5,'markSz',5);
% printas='-djpeg98';
printas=[];
figName='d:\hh\rapunzel';

% ----------------------------------------------------------------------
etslconst;
Etsl=[];
for g=1:size(fn,1)
  load([ddir fn(g,:)]);
  if g==1
    % if tOffset is nan set beginning of first recording as t=0
    if isnan(tOffset)
      tOffset=head.ds.fileInfo.recTime(1);
    end
  end
  bu.etsl(:,etslc.tsCol)=bu.etsl(:,etslc.tsCol)+(head.ds.fileInfo.recTime(1)-tOffset)*1000;
  Etsl=cat(1,Etsl,bu.etsl);
end



figure(1),hold on
orient portrait
plot(Etsl(:,etslc.tsCol)/1000,Etsl(:,etslc.durCol),'ro');
axis tight
set(gca,'xlim',[0 Etsl(end,etslc.tsCol)/1000]);
xlabel('time (s)')
ylabel('burst length (ms)')

% ----- printing -----
if ~isempty(printas), 
  print(printas,'-r300',figName); 
end
