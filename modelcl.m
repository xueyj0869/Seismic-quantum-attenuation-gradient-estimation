% This is a model test example
% 
%%%%%%%%%%%%%%%%%%%%%%%  model test    reservior thickness=15m
load('vmodel.mat');
figure;
imagesc(v')
xlabel('Depth(m)');
ylabel('Depth(m)');
title('Geological model with reservior thickness 15 m');
%%%%%%%%%%%%%%%%%%%%%  Import data and sampling frequency

% [seis,SegyTraceHeaders,SegyHeader]=ReadSegy('seismodel.sgy');
% seis=seis';
% Ts=SegyHeader.dt/1000000;
% fs=1/Ts;
%%%%%%%%%%%%%%%%%%%%%  Import data and sampling frequency
load('seis.mat');% seismic section for the geological model
fs=1000;
figure;
imagesc(seis');
xlabel('Trace');
ylabel('Time(ms)');
title('seismic section for the geological model');
out=quanattenver(seis,fs,0.018);
sigma1=5;  %%%%%高斯平滑
gauss_filter=fspecial('gaussian',15,sigma1);
Iout=imfilter(out,gauss_filter);
quanatt=(Iout-min(min(Iout)))./(max(max(Iout))-min(min(Iout)));

figure;
imagesc(quanatt');
xlabel('Trace');
ylabel('Time(ms)');
title('seismic attenuaton gradient estimation section');
load('mycolor.mat');
colormap(mycolor);colorbar
%%%%%%%%%%%%%%%%%% add noise  h smaller, f greater
% nseis=awgn(seis,10,'measured');
load('nseisSNR10.mat');
out=quanattenver(nseis,fs,0.01);   
sigma1=5;  %%%%%高斯平滑
gauss_filter=fspecial('gaussian',15,sigma1);   
Iout=imfilter(out,gauss_filter);
quanatt=(Iout-min(min(Iout)))./(max(max(Iout))-min(min(Iout)));

figure;
imagesc(quanatt');
xlabel('Trace');
ylabel('Time(ms)');
title('seismic attenuaton gradient estimation section');
load('mycolor.mat');
colormap(mycolor);colorbar
