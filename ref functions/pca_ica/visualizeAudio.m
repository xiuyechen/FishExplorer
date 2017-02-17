function visualizeAudio(y,name,Fs)
% Syntax:   visualizeAudio(y,name);
%           visualizeAudio(y,name,Fs);

% Constants
WINDOW_MS = 40;             % Spectrogram window width, in ms
Y_PAD     = 0.05;           % Ampitude (y-axis) padding

% Parse inputs
if ~exist('Fs','var') || isempty(Fs)
    % Default sampling rate
    Fs = 8000;
end

Ny    = numel(y);
t     = (1:Ny) / Fs;
axLim = [t(1), t(Ny), -(1 + Y_PAD), (1 + Y_PAD)];

% Compute spectrogram
window = WINDOW_MS * (Fs / 1000);
[S, F, T] = spectrogram(y,window,[],[],Fs,'yaxis'); 
P = 10 * log10(abs(S));

size(P)
size(F)
size(T)

% Plot waveform
subplot(2,1,1);
plot(t,y);
ax1 = gca();
axis(ax1,axLim);
axis(ax1,'manual');
hold(ax1,'on');
xlabel('Time (s)');
ylabel('Amplitude');
title('Waveform');

% Plot spectrogram
subplot(2,1,2);
pcolor(T,F,P);
shading flat;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram');

% Set figure name
fig = gcf();
set(fig,'Name',name);

% Play sound
sound(y);
