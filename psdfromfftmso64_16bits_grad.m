
function psdfromfftmso64_16bits(name,count,expNum,fs,save,folder)
mats = dir('mso*.mat');
% clearvars
% close all

% sample_pace = 1000;
gradient=1;
sgrad1=zeros(8000,1);sgrad2=sgrad1;
% auto file detecter:
% files = dir('NEWSUM2*.mat');
save_graphs = save;
sp = fs;
spt = sp*8;
%for mat = 1:size(files,1)

%filetitle = char(files(mat).name);
%load(filetitle);

load(name);
a=[1 3 5 7 9 11 13 15 17 19 21 23 25 ];b=[2 4 6 8 10 12 14 16 18 20 22 24 26];
n=1:8000;
Ptot=zeros(4001,1);
%plot(n,slopes(:,[a b]))
Gnum = 0;
for jj=1:(count/2)
    
disp(jj);
grad1=slopes(:,a(jj));
grad2=slopes(:,b(jj));
grad1=grad1-mean(grad1);
grad2=grad2-mean(grad2);
%sgrad1=sgrad1+grad1;
%sgrad2=sgrad2+grad2;
% figure
% plot(n,grad1,n,grad2)
% title('grad 1 and grad 2');
% xlabel('Mili Second');
% ylabel('Magnetic field');
fs=sp;
% f0=6.67;
% 
%  notchWidth=1e-3;
range1=max(grad1)-min(grad1);
range2=max(grad2)-min(grad2);
norm1=grad1/range1;norm2=grad2/range2;
% grad1=fnotch(grad1,fs,f0,notchWidth);
% grad2=fnotch(grad2,fs,f0,notchWidth);

% figure
% plot(n,norm1,n,norm2)
% title('normalized grad 1 and 2');
norm1=norm1/1;
norm2=norm2/1;
dif=(norm1-norm2)*min(range1,range2);
%dif=norm1*range1;
% plot(n,dif)

% fs=1000;
 %f0=4.66;
% % f0=50;
  %notchWidth=3e-2;
 %dif=fnotch(dif,fs,f0,notchWidth);
% f0=100;
% dif=fnotch(dif,fs,f0,notchWidth);
% f0=150;
% dif=fnotch(dif,fs,f0,notchWidth);
% f0=200;
% dif=fnotch(dif,fs,f0,notchWidth);
% f0=250;
% dif=fnotch(dif,fs,f0,notchWidth);

%  notchWidth=5e-4;
% f0=16.6;
% dif=fnotch(dif,fs,f0,notchWidth);
% f0=17.8;
% dif=fnotch(dif,fs,f0,notchWidth);
% figure
% plot(n,dif)
%[sp,vals]=spaps(n,dif,0.00000005);fdif=dif-0.95*vals';
%figure
%plot(n,dif,n,fdif);
fdif=dif;
% figure
% plot(n,dif);
% title('differences');
mf1=8.33e-7;
mf2=1/3.5e9;

mf=mf1/1;



Fs = sp;     % Sampling frequency                    
T = 1/Fs;      % Sampling period       
L = spt;      % Length of signal
t = (0:L-1)*T;        % Time vector
f = Fs*(0:(L/2))/L;

fig_singles = figure();
w= hann(spt); % window function
window= w.*sqrt(spt/sum(w.^2)); % normalize window for P= 1
%xw= x.*window';
fdif=fdif(1:spt).*window;
grad1=grad1(1:spt).*window;
grad2=grad2.*window;
% fdif=fdif(1:8000).*hann(8000);
% grad1=grad1(1:8000).*hann(8000);
% grad2=grad2.*hann(8000);

% grad1 = Notch_filter(grad1,fs,50,0.1);

A=fft(grad1);
%A=fft(grad1);
%B=fft(grad2);
PA2 = abs(A/L);
PA1 = PA2(1:L/2+1);
PA1(2:end-1) = 2*PA1(2:end-1);
loglog(f,PA1*mf/sqrt(1)/1)

filename = mats.name;
index = strfind(filename,'_');
filetitle = char(extractBetween(filename,1,index(2)-2));
filetitle = strcat(filetitle, num2str(expNum-1),char(extractBetween(filename,(index(2)),index(3)-1)));
filetitle = char(extractBetween(filetitle,1,index(3)-2));
filetitle = strcat(filetitle, num2str(Gnum),'_');
filetitle = strcat(filetitle,'_',' ch1 and ch3');
title(filetitle);


xlabel('Hz');
ylabel('Magnetic field / sqrHz');
grid on

hold on

% w= hann(8000); % window function
% window= w.*sqrt(8000/sum(w.^2)); % normalize window for P= 1
% %xw= x.*window';
% fdif=fdif(1:8000).*window;
% grad1=grad1(1:8000).*window;
% grad2=grad2.*window;
% fdif=fdif(1:8000).*hann(8000);
% grad1=grad1(1:8000).*hann(8000);
% grad2=grad2.*hann(8000);

% grad2 = Notch_filter(grad2,fs,50,0.1);

A=fft(grad2);
%A=fft(grad1);
%B=fft(grad2);
PA2 = abs(A/L);
PA1 = PA2(1:L/2+1);
PA1(2:end-1) = 2*PA1(2:end-1);
loglog(f,PA1*mf/sqrt(1)/1)
legend('magnetometer1','magnetometer2','1 cm gradiometer');


if save_graphs == 1
    filetitle = strcat(filetitle,'.pdf');
    saveas(fig_singles,fullfile(folder,filetitle));
end



fig_grad = figure();
% w= hann(8000); % window function
% window= w.*sqrt(8000/sum(w.^2)); % normalize window for P= 1
% %xw= x.*window';
% fdif=fdif(1:8000).*window;
% grad1=grad1(1:8000).*window;
% grad2=grad2.*window;
% fdif=fdif(1:8000).*hann(8000);
% grad1=grad1(1:8000).*hann(8000);
% grad2=grad2.*hann(8000);
% fdif = Notch_filter(fdif,fs,50,0.1);
A=fft(fdif);
%A=fft(grad1);
%B=fft(grad2);
PA2 = abs(A/L);
PA1 = PA2(1:L/2+1);
PA1(2:end-1) = 2*PA1(2:end-1);
loglog(f,PA1*mf/sqrt(1)/1)

% ----- making title -----------

filename = mats.name;
index = strfind(filename,'_');
filetitle = char(extractBetween(filename,1,index(2)-2));
filetitle = strcat(filetitle, num2str(expNum-1),char(extractBetween(filename,(index(2)),index(3)-1)));
filetitle = char(extractBetween(filetitle,1,index(3)-2));
filetitle = strcat(filetitle, num2str(Gnum),'_');
filetitle = strcat(filetitle, '_Gradient');
Gnum = Gnum + 1;
title(filetitle);
xlabel('Hz');
ylabel('grad mag field / sqrHz');
grid on
%PB2 = abs(B/L);


% PB1 = PB2(1:L/2+1);
% PB1(2:end-1) = 2*PB1(2:end-1);
% loglog(f,PB1*mf/sqrt(2)/1)
% PC2 = abs(C/L);
% PC1 = PC2(1:L/2+1);
% PC1(2:end-1) = 2*PC1(2:end-1);
% loglog(f,PC1*mf/sqrt(8))
% grid
legend('gradiometer','magnetometer2','1 cm gradiometer');
% xlabel('Frequency,(Hz)')
% ylabel('Magnetic Field ASD,(T/sqrt(Hz))')
% title('Magnetic Field Amplitude Spectral Density, 16 bits')
if save_graphs == 1
    filetitle = strcat(filetitle,'.pdf');
    saveas(fig_grad,fullfile(folder,filetitle));
end
end
%end


% srange1=max(sgrad1)-min(sgrad1);
% srange2=max(sgrad2)-min(sgrad2);
% snorm1=sgrad1/srange1; snorm2=sgrad2/srange2;
% sdif=(snorm1-snorm2)*min(srange1,srange2)/1;
% 
% 
% Fs = 1000;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L = 8000;             % Length of signal
% t = (0:L-1)*T;        % Time vector
% f = Fs*(0:(L/2))/L;
% figure
% sdif=sdif.*chebwin(8000);
% sgrad1=sgrad1.*chebwin(8000);
% sgrad2=sgrad2.*chebwin(8000);
% C=fft(sdif);
% A=fft(sgrad1);
% B=fft(sgrad2);
% PA2 = abs(A/L);
% PA1 = PA2(1:L/2+1);
% PA1(2:end-1) = 2*PA1(2:end-1);
% loglog(f,PA1*mf/sqrt(8)/1)
% hold on
% PB2 = abs(B/L);
% 
% 
% 
% 
% 
% PB1 = PB2(1:L/2+1);
% PB1(2:end-1) = 2*PB1(2:end-1);
% loglog(f,PB1*mf/sqrt(8)/1)
% PC2 = abs(C/L);
% PC1 = PC2(1:L/2+1);
% PC1(2:end-1) = 2*PC1(2:end-1);
% loglog(f,PC1*mf/sqrt(8))
% grid
% 
% 
% legend('magnetometer1','magnetometer2','30 cm gradiometer')
% xlabel('Frequency,(Hz)')
% ylabel('Magnetic Field ASD,(T/sqrt(Hz))')
% title('Magnetic Field Amplitude Spectral Density,TekMSO64 16 bits')
end

%fnotch function
function filRes = Notch_filter(slope,fs,f0,notchWidth)
    filRes = slope;
    for i=[1 2 3 4 5 7 9]
        filRes = fnotch(filRes,fs,(f0*i),notchWidth);
    end
end

        