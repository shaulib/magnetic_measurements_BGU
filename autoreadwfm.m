
function autoreadwfm()
clearvars
close all
% auto recognizing the number of files to work on:
numf = dir('*.wfm');
numfil = size(numf,1);
tic
filetitle='mso64';

% 20,24 creats mat, file. 20,20 create wfm files
names=zeros(20,24);
wfm=numf;
wfmcell=struct2cell(wfm);
namcell=wfmcell(1,1:numfil);
%creating mat files of samples:
 for i=1:numfil
    nam=cell2mat(namcell(1,i));
    names(i,:)=nam;
    savednam = [filetitle nam(1:20)];
    clear y;
    clear t;
    [y,t,info]=wfm2read(nam);
    save(savednam,'y')
    toc
 end
 toc
end


