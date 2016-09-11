%boiler plate for mtspecgramc
%define t and x
% t=.001:.001:30;
t=.004:.004:30;
f1=20; f2=50;
h = hamming(length(t)/2);
x= [h'.*sin(2*pi*f1*t(1:length(t)/2)) h'.*sin(2*pi*f2*t(length(t)/2+1:end))];

win_length = 1;
movingwin = [win_length win_length*.2];

dt = t(2)-t(1);
Fs = 1/dt;
T = win_length; %numel(t)/Fs;
W = 5; %.1;
TW=T*W;
ntapers = max(round(2*TW-1),1);
params.Fs = Fs;
params.tapers = [TW,ntapers];
params.pad = 0;
params.err = [2, 0.05];
params.fpass = [0 max(f1,f2)*2];
[Sxx,t,f,Serr]=mtspecgramc(x,movingwin, params);
imagesc(t, f, 10*log10(Sxx'))
% ylim([0 max(f1,f2)*2])
ylabel('Frequency [Hz]')
xlabel('Time [s]')
axis xy
colorbar

f_smear_each_side = TW/win_length
max_f_res = 1/win_length