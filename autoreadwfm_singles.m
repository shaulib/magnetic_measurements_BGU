clearvars
close all
% auto recognizing the number of files to work on:
numf = dir('*.wfm');
numfil = size(numf,1);
tic
filetitle='mso64';

names=zeros(20,20);
wfm=dir('*.wfm');
wfmcell=struct2cell(wfm);
namcell=wfmcell(1,1:numfil);
for i=1:numfil
    nam=cell2mat(namcell(1,i));
    names(i,:)=nam;
    nam
    savednam=[filetitle nam(1:20)]
    clear y;
    clear t;
    [y,t,info]=wfm2read(nam);
    save(savednam,'y')
    toc
end

% change files name to .mat
files = dir('mso*.wfm');
for id = 1:length(files) 
    [~, f,ext] = fileparts(files(id).name);
    rename = strcat(f,'','.mat'); 
    movefile(files(id).name, rename); 
end

toc
