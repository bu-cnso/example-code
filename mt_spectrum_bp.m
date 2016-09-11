%boiler plate for mtspectrumc
%define t and x
t=.001:.001:5;
f1=20; f2=50;
x=sin(2*pi*f1*t)+sin(2*pi*f2*t);
dt = t(2)-t(1);
Fs = 1/dt;
T = numel(t)/Fs;
W = 5;
TW=T*W;
ntapers = round(2*TW-1);
params.Fs = Fs;
params.tapers = [TW,ntapers];
params.pad = 0;
params.err = [2, 0.05];
[Sxx,f,Serr]=mtspectrumc(x,params);

% plot(f, 10*log10(Sxx), 'r') %Plot power vs frequency.
% hold on
% plot(f, 10*log10(Serr), 'b')
% hold off

subplot(211)
plot(t,x)
subplot(212)
shadedErrorBar(f,10*log10(Sxx),10*log10(Serr), 'b')

xlabel('Frequency [Hz]')
ylabel('Power [ mV^2/Hz]')

f_smear_each_side = W % = TW/T = W