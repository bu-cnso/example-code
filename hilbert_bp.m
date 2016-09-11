%boiler plate for hilbert transform
%define t and x
t=.001:.001:1;
% f1=20; f2=30;
% x=sin(2*pi*f1*t)+sin(2*pi*f2*t);
f = 25;
x=sin(2*pi*f*t);

dt = t(2)-t(1);
Fs = 1/dt;
T = numel(t)/Fs;

fpass = [20,30];
Wn = [fpass(1)*2/Fs,fpass(2)*2/Fs];               %Define the low frequency window of interest.
[B,A] = butter(2,Wn,'bandpass');  %Apply the filter to isolate the band.
filtered_signal = filtfilt(B,A,x);

phi = angle(hilbert(filtered_signal));    %Phase
amp = abs(hilbert(filtered_signal));      %Amplitude envelope

ax(1) = subplot(2,1,1);
plot(t, x,'LineWidth', 2);
hold on
ax(2) = axes('Position',ax(1).Position,...
    'YAxisLocation','right',...
    'Color','none');
plot(ax(2),t, phi, 'g','LineWidth', 2);
set(ax(2),'Position',ax(1).Position,...
    'YAxisLocation','right',...
    'Color','none');
hold off
axis tight
xlabel('Time [s]');  title('Signal and Phase');

ax(3) = subplot(2,1,2);
plot(t, x,'LineWidth', 2)
hold on
plot(t, amp, 'r','LineWidth', 2)
hold off
axis tight
xlabel('Time [s]');  title('Signal and amplitude envelope')

linkaxes(ax,'x')

%% angles

meanPhi = angle(mean(exp(1i*phi)));

figure
subplot(121)
rPlot = rose(phi);
hold on
% xR = get(rPlot,'Xdata');
% yR = get(rPlot,'Ydata');
% patch(xR,yR,'b');
compass(exp(1i*meanPhi) * .65*sqrt(max(get(gca,'xlim'))^2 +max(get(gca,'ylim'))^2)  , 'r')
title('Phase Roseplot');

subplot(122)
hold on
hist(phi,20)
ylim = get(gca,'ylim');
plot([meanPhi meanPhi], ylim, 'r', 'Linewidth',2)
hold off
xlim([-pi,pi])
title('Phase Histogram');
