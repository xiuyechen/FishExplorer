%%
U = unique(gIX);
numU = length(U);
C = zeros(numU,size(M,2));
D = zeros(numU,1);
for i=1:numU,
    IX = find(gIX == U(i));
    if length(IX)==1,
        C(i,:) = M(IX,:);
        D(i) = 1;
    else
        M_s = M(IX,:);
        [~,C1,~,D1] = kmeans(M_s,1,'distance','correlation');
        C(i,:) = C1;
        D(i) = mean(D1);
    end
end
%%

figure
Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = 1000;                     % Length of signal
t = (0:L-1)*T;                % Time vector
% Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
x = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t); 
y = x + 2*randn(size(t));     % Sinusoids plus noise
plot(Fs*t(1:50),y(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('time (milliseconds)')

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

%%
Fs = 1.97;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = size(C,2);                     % Length of signal

H = zeros(1,size(C,1));
for i =1:size(C,1),
y = C(i,:);

NFFT = 2^nextpow2(L); 
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
Amp = 2*abs(Y(1:NFFT/2+1));

thres_peakheight = prctile(Amp,99);
[pks,locs] = findpeaks(Amp,'MinPeakDistance',5,'MinPeakHeight',thres_peakheight);
prd = 1./f(locs);
[~,I] = min(abs(locs-80));
H(i) = pks(I);
end

figure;plot(H)
% findpeaks(x, 'MINPEAKDISTANCE', dist);
% figure; hold on;
% plot(f,Amp);plot(f(locs),Amp(locs),'ro')
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')

% 1./f(locs)
%%
% Fs = 1;                    % Sampling frequency
% T = 1/Fs;                     % Sample time
% L = 1000;                     % Length of signal
% t = (0:L-1)*T;                % Time vector
%  y = zeros(1,1000);y(1:10:1000) = 1;
% 
% NFFT = 2^nextpow2(L); 
% Y = fft(y,NFFT)/L;
% f = Fs*linspace(0,1,NFFT);
% Amp = abs(Y(1:NFFT));
% 
% thres_peakheight = prctile(Amp,99);
% [pks,locs] = findpeaks(Amp,'MinPeakDistance',5,'MinPeakHeight',thres_peakheight);
% % findpeaks(x, 'MINPEAKDISTANCE', dist);
% figure; hold on;
% plot(f,Amp);plot(f(locs),Amp(locs),'ro')
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')

%%
Fs2 = 1; %?
L2 = size(Amp,2);
NFFT2 = 2^nextpow2(L2); 
Y2 = fft(Amp,NFFT2)/L2;
f2 = Fs2/2*linspace(0,1,NFFT2/2+1);
Amp2 = 2*abs(Y2(1:NFFT2/2+1));
figure;
plot(f2,Amp2);
