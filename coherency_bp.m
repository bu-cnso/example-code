%boiler plate for coherencyc
clear

% dbstop if error
%define t and x
noiseAmplitude = 1;
t=.001:.001:5;
f1=20; f2=50;

num_trials = 10;
shift = rand(1, num_trials);
% shift = repmat(1/f2/10,num_trials,1);
% shift = repmat(0,num_trials,1);
for trial = 1:num_trials
    x1(trial,:)=sin(2*pi*f1*t)+sin(2*pi*f2*t)+noiseAmplitude * rand(1, length(t));
    x2(trial,:)=sin(2*pi*f1*t+shift(trial))+sin(2*pi*f2*t-shift(trial))+noiseAmplitude * rand(1, length(t));
end
tic
dt = t(2)-t(1);
Fs = 1/dt;
T = numel(t)/Fs;
W = 5;
TW=T*W;
ntapers = 2*TW-1;
params.Fs = Fs;
params.tapers = [TW,ntapers];
params.pad = 0;
params.err = [2, 0.05];
if num_trials>1
    params.trialave = 1;
else
    params.trialave = 0;
end
params.fpass = [0 100];

x1_0mean = bsxfun(@minus,x1,mean(x1,2))';
x2_0mean = bsxfun(@minus,x2,mean(x2,2))';
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(x1_0mean, x2_0mean, params);

subplot(311)
plot(f, real(S1), 'r') %Plot power vs frequency.
subplot(312)
plot(f, real(S2), 'b')
subplot(313)
Cerr2 = [Cerr(2,:);Cerr(1,:)];
shadedErrorBar(f,(C),Cerr2, 'b')
% plot(f, C)

xlabel('Frequency [Hz]')
ylabel('Coherence')

f_smear_each_side = W % = TW/T = W
if num_trials>1
    toc/num_trials
else
    toc
end
%% timing
%.32 per single
%.75 per single with 10
%1.5 per single with 20