clc, clearvars, close all;
fs = 10*2000; %sampling frequency
t = 0:1/fs:4;
% 1
f1 = 500;
f2 = 1000;
f3 = 1500;
f4 = 2000;
x_t = cos(2*pi*f1*t) + cos(2*pi*f2*t) + cos(2*pi*f3*t) + cos(2*pi*f4*t);
% 2
audiowrite('xt.wav', x_t, fs);
% 3
figure;
plot(t, x_t);
title('signal x(t)');
xlabel('t (seconds)');
ylabel('x(t)');
% 4
E_t = sum(abs(x_t).^2) * (1/fs);
% 5
X_f = fft(x_t);
% 6
N = length(X_f); % Number of samples
f = (-fs/2):(fs/N):(fs/2-fs/N);
X_f_shifted = fftshift(X_f);
figure;
plot(f,abs(X_f_shifted));
title('magnitude of X(f)');
xlabel('frequency (HZ)');
ylabel('|X(f)|');
% 7
E_freq = (sum(abs(X_f).^2) / N) * (1/fs);
disp(['Time-domain energy E_t: ', num2str(E_t)]);
disp(['Frequency-domain energy E_freq: ', num2str(E_freq)]);
% Check if energies match
if abs(E_t - E_freq) < 1e-10 % Tolerance for numerical precision
    disp('Parseval''s theorem verified');
else
    disp('Parseval''s theorem not verified');
end
% 8
filter_order = 20;
f_cutoff = 1250;
wn = f_cutoff / (fs/2); % normalized (0 to 1 scale)
[low_b, low_a] = butter(filter_order , wn , "low");
% 9
figure;
freqz(low_b, low_a, fs, fs);
% 10
y1_t = filter(low_b, low_a, x_t);
% 11
audiowrite('y1_t.wav', y1_t, fs);
% 12
figure;
plot(t,y1_t);
title('Signal y1(t)');
xlabel('time (seconds)');
ylabel('y1(t)');
% 13
Ey1_t = sum(abs(y1_t).^2) * (1/fs);
% 14
Y1_f = fft(y1_t);
% 15
N2 = length(Y1_f); % Number of samples
f = (-fs/2):(fs/N2):(fs/2-fs/N2);
Y1_f_shifted = fftshift(Y1_f);
figure;
plot(f,abs(Y1_f_shifted));
title('magnitude of Y1(f)');
xlabel('frequency (HZ)');
ylabel('|Y1(f)|');
% 16
Ey1_freq = (sum(abs(Y1_f).^2) / N2) * (1/fs);
disp(['Time-domain energy Ey1_t: ', num2str(Ey1_t)]);
disp(['Frequency-domain energy Ey1_freq: ', num2str(Ey1_freq)]);
% Check if energies match
if abs(Ey1_t - Ey1_freq) < 1e-10 % Tolerance for numerical precision
    disp('Parseval''s theorem verified');
else
    disp('Parseval''s theorem not verified');
end
% 17
f_cutoff_hp = 1250; 
wn_hp = f_cutoff_hp / (fs/2); % Normalized frequency
[high_b, high_a] = butter(filter_order, wn_hp, "high");
sos_high = tf2sos(high_b, high_a);
% 18
figure;
freqz(high_b, high_a, fs, fs);
title('Butterworth High-Pass Filter: Magnitude and Phase Response');

% 19
% Apply the signal x(t) to the Butterworth HPF
y2_t = sosfilt(sos_high,x_t);

% 20
audiowrite('y2_t.wav', y2_t, fs);

% 21
% Plot the signal y2_t versus time t
figure;
plot(t, y2_t);
title('Signal y2(t) (High-Pass Filtered)');
xlabel('Time (seconds)');
ylabel('y2(t)');

% 22
Ey2_t = sum(abs(y2_t).^2) * (1/fs);
disp(['Time-domain energy Ey2_t: ', num2str(Ey2_t)]);

% 23
Y2_f = fft(y2_t);

% 24
% Plot the magnitude of Y2(f) in the frequency range -fs/2 to fs/2
N3 = length(Y2_f); % Number of samples
f = (-fs/2):(fs/N3):(fs/2-fs/N3);
Y2_f_shifted = fftshift(Y2_f);
figure;
plot(f, abs(Y2_f_shifted));
title('Magnitude of Y2(f)');
xlabel('Frequency (Hz)');
ylabel('|Y2(f)|');

% 25
Ey2_freq = (sum(abs(Y2_f).^2) / N3) * (1/fs);
disp(['Frequency-domain energy Ey2_freq: ', num2str(Ey2_freq)]);

% Verify Parseval's theorem
if abs(Ey2_t - Ey2_freq) < 1e-10 % Tolerance for numerical precision
    disp('Parseval''s theorem verified for y2(t)');
else
    disp('Parseval''s theorem not verified for y2(t)');
end

