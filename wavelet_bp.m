%% wavelet_bp
% boiler plate for wavelet spectrum

fs = 500; %hz
T = 30; % s
t=1/fs:1/fs:T;
f1=20; f2=50;
h = hamming(length(t)/2);
x= [h'.*sin(2*pi*f1*t(1:length(t)/2)) h'.*sin(2*pi*f2*t(length(t)/2+1:end))]';
N = length(x);
x = x+normrnd(0,.5,N,1);

%% from ben
% freqs = 1:100;
% nFreqs = length(freqs);
% nCycles = linspace(3, 25, nFreqs);
% 
% bands = [1 4; 4 8; 8 30; 30 100];
% noBands = size(bands, 1);
% 
% cycle_lengths = nCycles.*(sampling_freq./freqs);
% 
% wavelets = dftfilt3(freqs, nCycles, sampling_freq, 'winsize', max(sampling_freq, max(cycle_lengths)));
% segment_length = sampling_freq;
% [Spec, Spec_norm, Spec_pct, Spec_norm_pct] = deal(nan(length(PD_dec), nFreqs, 2));
% data = PD_dec(:, ch);
% data_reflected = [flipud(data(1:segment_length)); data; flipud(data((end - segment_length + 1):end))];
% for f = 1:nFreqs
%   conv_prod = conv(data_reflected, wavelets(f,:), 'same');
%   Spec(:, f, ch) = conv_prod((segment_length + 1):(end - segment_length));
% end
% 
% % Baseline normalize.
% baseline_mean = mean(abs(Spec(t <= basetime, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
% Spec_pct(:, :, ch) = 100*abs(Spec(:, :, ch))./(ones(size(Spec(:, :, ch)))*diag(baseline_mean)) - 100;
% 
% % Normalize by total power.
% Spec_norm(:, :, ch) = abs(Spec(:, :, ch))./repmat(sqrt(sum(abs(Spec(:, :, ch)).^2, 2)), 1, nFreqs);
% 
% % Baseline normalize percent of total power.
% baseline_mean = mean(abs(Spec_norm(t <= basetime, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
% Spec_norm_pct(:, :, ch) = 100*abs(Spec_norm(:, :, ch))./(ones(size(Spec_norm(:, :, ch)))*diag(baseline_mean)) - 100;

%% book
% s = nCycles/(2*pi*f); %std of gaussian taper
% gaussian_win = exp(-time.^2./(2*s^2));

nFreqs      =  30;  % number of frequency bands
lowestFrequency  =   2;  % in Hz
highestFrequency = 100;  % in Hz
lowestCycles = 3;
highestCycles = 25;

cyclesForFreq = linspace(lowestCycles, highestCycles, nFreqs);
% cyclesForFreq = logspace(log10(lowestCycles), log10(highestCycles), nFreqs);
waveletFrequencies = linspace(lowestFrequency, highestFrequency, nFreqs);
% waveletFrequencies = logspace(log10(lowestFrequency), log10(highestFrequency), nFreqs);

waveletTimeLengthSec = 2*cyclesForFreq(1)/lowestFrequency;
waveletTime = (-waveletTimeLengthSec/2 : 1/fs : waveletTimeLengthSec/2)'; % time span of wavelets
waveletLength = length(waveletTime);
waveletHalfLength = ceil(waveletLength/2);

nConv = waveletLength + N - 1;
nConvPow2 = pow2(nextpow2(nConv));
fftData = fft(x, nConvPow2);

s = cyclesForFreq./(2*pi*waveletFrequencies);
A = sqrt(1./(s*sqrt(pi)));

% % initialize power for time x freq
% dataSpectrogram = zeros(N, nFreqs);
% 
% % init power for freq
% dataSpectrum = zeros(nConvPow2, nFreqs);
% 
% % initialize wavelet mats
% waveletFamily = zeros(waveletLength, nFreqs);
% waveletPower = zeros(nConvPow2, nFreqs);

% %Loop through frequencies and make a family of wavelets.
% tic
% for iFreq = 1:nFreqs
% %     Explanation:
% %     % create a sine wave at this frequency
% %     sinewave = exp(2*1i*pi * frequencies(iFreq).*waveletTime); % the "1i" makes it a complex wavelet
% %     
% %     % create a Gaussian window
% %     gaus_win = exp(-waveletTime.^2./(2*(cyclesForFreq(iFreq)/(2*pi*frequencies(iFreq)))^2));
% %     
% %     % create wavelet via element-by-element multiplication of the sinewave and gaussian window
% %     waveletFamily(iFreq,:) = sinewave.*gaus_win;
%     
%   thisWavelet = A(iFreq)*exp(2*1i*pi*waveletFrequencies(iFreq).*waveletTime - waveletTime.^2./(2*s(iFreq)^2));
%   waveletFamily(:, iFreq) = thisWavelet;
%   
%   thisFftWavelet = fft(thisWavelet, nConvPow2);
%   waveletPower(:, iFreq) = thisFftWavelet .* conj(thisFftWavelet);
%   
%   thisFftWaveletData = thisFftWavelet.*fftData;
%   dataSpectrum(:, iFreq) = thisFftWaveletData .* conj(thisFft WaveletData);
%   
%   dataWavConv = ifft(thisFftWaveletData);
%   dataWavConv = dataWavConv(1:nConv);
%   dataWavConv = dataWavConv(waveletHalfLength+1:end-waveletHalfLength+2);
%   
%   dataSpectrogram(:,iFreq) = abs(dataWavConv).^2; %time x freq
% end
% toc

%vectorized version of above
% tic
waveletFamily = bsxfun(@times, A, exp( bsxfun(@times, 2*1i*pi*waveletFrequencies, waveletTime)...
  - bsxfun(@rdivide, waveletTime.^2, (2*s.^2)) ));
fftWavelet = fft(waveletFamily, nConvPow2);
waveletPower = fftWavelet .* conj(fftWavelet);
fftWaveletData = bsxfun(@times, fftWavelet, fftData);
dataSpectrum = fftWaveletData .* conj(fftWaveletData);
dataWavConv = ifft(fftWaveletData);
dataWavConv = dataWavConv(1:nConv, :);
dataWavConv = dataWavConv(waveletHalfLength+1:end-waveletHalfLength+2, :);
dataSpectrogram = abs(dataWavConv).^2;
% toc

% fttWD = bsxfun(@times, fftData, fftWavelet);
% sxx = fttWD.* conj(fttWD);

df=fs/nConvPow2;
wavFftFreq=0:df:fs-df;
lastF = floor(nConvPow2/2) - 1;
wavFftFreq = wavFftFreq(1:lastF);
waveletPower = waveletPower(1:lastF,:);
dataSpectrum = dataSpectrum(1:lastF,:);

%% Wavelet Power Plot
figure
s(1) = subplot(211);
plot(wavFftFreq, waveletPower)
xlim([0 waveletFrequencies(end)+10])
title('Wavelet Power 2d')
set(gca, 'xtick',waveletFrequencies)
set(gca, 'xticklabel',round(waveletFrequencies))
xlabel('Center Freq [Hz]')
ylabel('Power')
box off

s(2) = subplot(212);
imagesc(wavFftFreq,waveletFrequencies,waveletPower')
axis xy;
xlim([0 waveletFrequencies(end)+10])
set(gca, 'xtick',waveletFrequencies)
set(gca, 'xticklabel',round(waveletFrequencies))
title('Wavelet Power 3d')
xlabel('Center Freq [Hz]')
ylabel('Freq span [Hz]')

linkaxes(s,'x')

%% Data Plot
figure
s(1) = subplot(311);
plot(t,x)
title('Raw Data')
xlabel('Time [s]')
ylabel('Amplitude')
box off

s(2) = subplot(312);
imagesc(t,waveletFrequencies,dataSpectrogram')
axis xy;
title('Data Power Spectrogram')
xlabel('Time [s]')
ylabel('Freq [Hz]')

subplot(313)
plot(wavFftFreq, max(dataSpectrum,[],2))
title('Data Power Spectrum')
xlabel('Freq [Hz]')
ylabel('Power')

linkaxes(s,'x')