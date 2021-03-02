% clearvars
% close all
function cleverparse3125MSs16bits(fs)
tic
sav=1;
mats = dir('mso*.mat');
lastName = mats(size(mats,1)).name;
num_dex = strfind(lastName,'_');
exp_num = str2num(cell2mat(extractBetween(lastName,num_dex(1)+1,num_dex(2)-1)));
sp = fs;
spt = sp*8;

for exp = 0:exp_num

filename = mats.name;
index = strfind(filename,'_');
if exp > 9
    filetitle = char(extractBetween(filename,1,index(2)-3));
    filetitle = strcat(filetitle, num2str(exp),'_');
else
filetitle = char(extractBetween(filename,1,index(2)-2));
filetitle = strcat(filetitle, num2str(exp),'_');
end

% filetitle='mso6420200916_001_';
ln1=length(filetitle);
filehead=strcat(filetitle,'*');
kik=length(filetitle)+24;

files=dir(filehead);
kjk=length(files);
names=zeros(kjk,kik);
filescell=struct2cell(files);
ncell=filescell(1,:);
lfiles=length(files);
Sr=31.25e6;nomf=sp;



le=250000000;nv=50;numv=floor(le/nv);tlast=floor(le/Sr);
vary=zeros(1,numv);numw=tlast*nomf;stwf=zeros(1,numw);
lvar=Sr/nomf/nv;
slopes=zeros(numw,lfiles);
load('filtercoeff.mat')
load('bandp')
jstop=1;
  for j11=1:lfiles
%for j11=1:1
     nam=cell2mat(ncell(1,j11));
    
    nam
    st=ln1+1;fin=ln1+7;
   loadnam(j11,:)=[filetitle nam(st:fin)]
    
    clear y;
   
    load (loadnam(j11,:),'y');
    
while jstop==1

for i=1:numv
    low=(i-1)*nv+1;
    hi=i*nv;
    yy=y(low:hi);
    %tt(i)=t(low);
      vary(i)=var(yy);
end




for i=1:numw
  for k=1:lvar;
    if vary(lvar*(i-1)+k)>0.02
        stwf(i)=(i-1)*lvar*nv+nv*(k-1)+1;break
  end
    end
end
jstop=0;
end

lnspec=1;err=zeros(1,numw);db=err;longbigdb=zeros(1,numw*lnspec);
 bigerr1=zeros(lnspec,numw);  bigerr=zeros(lnspec,numw);bigerr0=bigerr1;bigdb=bigerr;
 bigdb0=bigerr0;bigdb1=bigerr1;  
  WS=71;
 bb=(1/WS)*ones(1,WS);
 ka=1;
mfsig11=zeros(1,spt);az=zeros(1,100);

b=aa;ka=1;

for i=1:nomf*tlast
 
    k1=stwf(i)+2001;k2=k1+20000;  ccc=1:2:40001;cc=1:40000;
   sig0=(y(k1:k2));sig1=interp1(ccc,sig0,cc);
   fsig1=filter(bb,ka,sig1);fsig1a=filter(bb,ka,fsig1);
%    fsig11=fsig1a;
L=length(sig1);
Bsig1=reshape(fsig1a,[5 L/5]);
mBsig1=mean(Bsig1);

%mfsig11=filter(b,ka,mBsig1);
mfsig11b=filter(d,mBsig1);
mfsig11=mfsig11b(2001:end);
%mfsig11=filter(bb,ka,mBfsig1);

% k=0
% 
% for j=1:5:L
%     k=k+1;
%     mfsig11(k)=mean(fsig11(j:j+4));
%    
%    
% end

    meana=mean(mfsig11(1:200));
    meanb=mean(mfsig11(end-200:end));
    mdif=meanb-meana;xlen=length(mfsig11);
    xc=1:xlen;xlen1=xlen-1;
    mfsig11=mfsig11-xc*mdif/xlen1+mdif/xlen1-meana;
if i==1
    figure
    plot(mfsig11)
end
   if i==4000
       figure
       plot(mfsig11)
end
ma=max(mfsig11);
mi=min(mfsig11);
msig=(mfsig11-mi)/(ma-mi)-0.5;


LL=length(msig);
ff1=msig(1:LL-1);
ff2=msig(2:LL);
ff=ff1.*ff2;

a=ff<0;az=find(a,100);azz=zeros(1,length(az));
ax1=1:length(az);
for j=1:length(az)
    ind=az(j);
    azz(j)=az(j)+abs(msig(ind))/(abs(msig(ind))+abs(msig(ind+1)));
   
end
format long;
AX1=[ones(length(ax1),1) ax1'];
   
b1=AX1\azz';
db(i)=b1(2);
if i==2500;toc
end


end
slopes(:,j11)=db;
end

toc
figure
plot(slopes)
if sav==1 ;save(['NEWSUM2' filetitle],'slopes');
end
end
end
    
    