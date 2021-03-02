  
% calling the 3 procedures read, parse and fft
sample_rate = 1000;
save_graphs = 1;
folder = 'D:\shauli\graph';
pairs_grad = 1;

% autoreadwfm_singles();
autoreadwfm();
cleverparse3125MSs16bits(sample_rate);
 
% auto file detecter:
files = dir('NEWSUM2*.mat');


for mat = 1:size(files,1)
filetitle = char(files(mat).name);
curfile = load(files(mat).name);
graphNum = size(curfile.slopes,2);
  if pairs_grad == 1
    psdfromfftmso64_16bits_grad(filetitle,graphNum,mat,sample_rate,save_graphs,folder);
    else
    psdfromfftmso64_16bits(filetitle,graphNum,mat,sample_rate,save_graphs,folder);
  end
end