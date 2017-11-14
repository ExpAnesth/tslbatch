% plots time course of firing rate for a given set of data

ddir='F:\ll\rawData\2009_01_08_set3\';
fIndex=0:16;

numOfFiles=length(fIndex);
fn=[];
% from the vector of integers make a vector of strings with width 4, namely
% the file indices (zero-padding)
for i=1:numOfFiles
  fn(i,:)=['2009_01_08_' sprintf('%.4i',fIndex(i))  '_IN0_SPX_res.mat']; 
end;

% graphics options
labelscale('fontSz',12,'scaleFac',1,'lineW',.5,'markSz',10);
printas=[];
% printas='-djpeg98';
figName='F:\ll\figures_report\verlauf';

binW=5000;

% ----------------------------------------------------------------------

Tsl=[];
for g=1:size(fn,1)
  load([ddir fn(g,:)]);
  if g==1
    tOffset=head.ds.fileInfo.recTime(1);
  end
  Tsl=cat(1,Tsl,evt.tsl{1}+(head.ds.fileInfo.recTime(1)-tOffset)*1000);
end
bin=0:binW:Tsl(end);

n=histc(Tsl,bin);


figure(1), clf
orient portrait
rexy('ax',gca,'xfac',1.0,'yfac',.5)
bar(bin/1000,n/binW*1000,1.0,'k')
niceyuax;
xlabel('time (s)')
ylabel('firing rate (Hz)')

% ----- printing -----
if ~isempty(printas), 
  print(printas,'-r500',figName); 
end
