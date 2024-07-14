% Discrete Fourier transform using fft function
% 08/12/2021 R. Liupekevicius
clear all; close all;
clc;

% the following example is found in mathworks doc fft

Fs = 1500;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector

% % test this function utility
% n = 2^nextpow2(L);

% Form a signal containing a 50 Hz sinusoid of amplitude 0.7 and a 120 Hz
% sinusoid of amplitude 1.
s50  = 0.7*cos(2*pi*50*t);
s120 = sin(2*pi*120*t);
S    = s50+s120;

% Corrupt the signal with zero-mean white noise with a variance of 4.
X = S + 2*randn(size(t));

%Plot the noisy signal in the time domain. It is difficult to identify
% the frequency components by looking at the signal X(t).
figure(1)
plot(1000*t(1:50),X(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')

%Plot full noisy signal
figure(2)
plot(1000*t,X)
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')

%% FFT

%Compute the Fourier transform of the signal.
Y = fft(X);

%Compute the two-sided spectrum P2. Then compute the single-sided spectrum
% P1 based on P2 and the even-valued signal length L.
P2 = abs(Y/L);   %scaling that takes place for all freq.
P1 = P2(1:L/2+1);%select only positive frequencies
P1(2:end-1) = 2*P1(2:end-1);%scalling for all freqs but f=0.

%Define the frequency domain f and plot the single-sided amplitude
% spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected,
% because of the added noise. On average, longer signals produce better
% frequency approximations.
f = Fs*(0:(L/2))/L;

figure(3)
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t), P1')
xlabel('f (Hz)')
ylabel('|P1(f)|')


%Plot real and imaginary parts of fft

Y = Y(1:L/2+1);
Y = 2*Y/L;
figure(5)
subplot(2,1,1)
plot(f,real(Y));
xlabel('f (Hz)')
ylabel('real 2Y/L')
subplot(2,1,2)
plot(f,imag(Y));
xlabel('f (Hz)')
ylabel('imag 2Y/L')
title('Single-Sided Amplitude Spectrum of X(t), 2Y/L')

%% Now, take the Fourier transform of the original, uncorrupted signal and
% retrieve the exact amplitudes, 0.7 and 1.0.

Yp = fft(S); %pure signal
P3 = abs(Yp/L); %scaling
P4 = P3(1:L/2+1); %taking only positive frequencies
P4(2:end-1) = 2*P4(2:end-1); %only multiply by 2 frequencies dif than zero

figure(4)
plot(f,P4) 
title('Single-Sided Amplitude Spectrum of S(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


% Yp = Yp(1:L/2+1);
% Yp = 2*Yp/L;
% figure(6)
% subplot(2,1,1)
% plot(f,real(Yp));
% xlabel('f (Hz)')
% ylabel('real 2Y/L')
% subplot(2,1,2)
% plot(f,imag(Yp));
% xlabel('f (Hz)')
% ylabel('imag 2Y/L')


fp = (0:L-1)/L*Fs - Fs/2;
Yp = 2*Yp/L;

figure(6)
subplot(2,1,1)
plot(fp,real(fftshift(Yp)));
xlabel('f (Hz)')
ylabel('real 2Y/L')
title('DFT');
subplot(2,1,2)
plot(fp,imag(fftshift(Yp)));
xlabel('f (Hz)')
ylabel('imag 2Y/L')




%% Using the hilbert transform on the time signal, then fft

%build frequency axis for fft
fh = (0:L-1)/L*Fs - Fs/2;


Sh=s50;
% Sh=hilbert(Sh);
Yh=fft(Sh);
Yh=2*Yh/L;

figure(7)
subplot(2,1,1)
plot(fh,real(Yh));
xlabel('f (Hz)')
ylabel('real 2Y/L')
subplot(2,1,2)
plot(fh,imag(Yh));
xlabel('f (Hz)')
ylabel('imag 2Y/L')
title("Y=fft(input time signal)");


%Explore fftshift function
Yh=fftshift(Yh);


figure(8)
subplot(2,1,1)
plot(fh,real(Yh));
xlabel('f (Hz)')
ylabel('real 2Y/L')
subplot(2,1,2)
plot(fh,imag(Yh));
xlabel('f (Hz)')
ylabel('imag 2Y/L')
title("Y=fft(input time signal)");

