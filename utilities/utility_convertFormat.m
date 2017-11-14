% saves analyzed data in slightly different (updated) format and slightly
% different (updated) file name

job=1;
% generate list of file names including full path of the files to be modified
switch job
  case 1
    cd d:\public\dr\Spike-Detection
    !dir *_tsl.mat /a /s /b > curFList.txt
end

% read that file, each file listed being assigned to a cell element
fili=textread('curFList.txt','%s','delimiter','\n');
disp(strvcat(fili));
bn=questdlg('contents of *.mat files as listed in command window will be saved as new variables and in a new file - sure you want to do this?','A little intimidation');

if strcmpi(bn,'yes')
  for ii=1:length(fili)
    switch job
      case 1
        % initialize new vars
        bu.etsl=[];
        bu.silentEtsl=[];
        evt.tsl=[];
        head=[];
        % new file name
        newFn=strrep(fili{ii},'_tsl.mat','_spx_res.mat');
        % load..
        load(fili{ii});
        % assign (spx: bu fields must be empty)
        evt.tsl=tsl;
        bu.etsl=[];
        bu.silentEtsl=[];
        head.ds.fileInfo=abfi;
        % save
        save(newFn,'evt','bu','head');
    end
  end
end

        
