% Discrete Fourier transform using fft function CHECK transfer function
% 08/12/2021 R. Liupekevicius
clear all; close all;
clc;


% definitions
L = 2048; %signal length
Fs= 2048; %samplig frequency
Ts = 1/Fs; %period of time steps
t = (0:(L-1))*Ts; %time vector;

% signal definition
mic1 = 0.7*sin(2*pi*50*t);
mic2 = sin(2*pi*50*t+ pi/3);

%fft both signals
Y1 =fft(mic1);
Y2 =fft(mic2);

% % fft both signals
% Y1 =fft(hilbert(mic1));
% Y2 =fft(hilbert(mic2));

% build frequency vector positive frequencies
fp = (0:L/2-1)/L*Fs;

% build frequency vector positive &n negative freqs
f  = (0:L-1)/L*Fs -Fs/2;

%% plot discrete fourier transform (DFT) on positive freq axis
figure(1);
hold on;
plot(fp,abs(2*Y1(1:L/2)/L),"linewidth", 1.5,"LineStyle","--", "Color","r");
plot(fp,abs(2*Y2(1:L/2)/L), "Color", "b");
legend();
xlabel('frequency[Hz]');
grid on;

%% plot DFT real and imag parts
figure(2)
title('real and imag parts of two signals with phase diff')
subplot(2,1,1)
hold on;
plot(f,real(2/L*fftshift(Y1)),"linewidth", 1.5,"LineStyle","--", "Color","r");
plot(f,real(2/L*fftshift(Y2)), "Color", "b");
hold off;
xlabel('frequency[Hz]');
legend();

subplot(2,1,2)
hold on;
plot(f,imag(2/L*fftshift(Y1)),"linewidth", 1.5,"LineStyle","--", "Color","r");
plot(f,imag(2/L*fftshift(Y2)), "Color", "b");
hold off;
xlabel('frequency[Hz]');
legend();


%% plot DFT abs and phase
figure(3)
title('abs and phase of two signals with phase diff')
subplot(2,1,1)
hold on;
plot(f,abs(2/L*fftshift(Y1)),"linewidth", 1.5,"LineStyle","--", "Color","r");
plot(f,abs(2/L*fftshift(Y2)), "Color", "b");
hold off;
xlabel('frequency[Hz]');
legend();

subplot(2,1,2)
hold on;
plot(f,180/pi*(angle(2/L*fftshift(Y1))),...
    "linewidth", 1.5,"LineStyle","--", "Color","r");
plot(f,180/pi*(angle(2/L*fftshift(Y2))), ...
    "Color", "b");
hold off;
xlabel('frequency[Hz]');
legend();


%% compute transfer function between mic1 and mic2


H12 = Y1./Y2;

figure(4)
subplot(2,1,1)
hold on;
plot(f,real(fftshift(H12)),"linewidth", 1.5,"LineStyle","--", "Color","r");
hold off;
xlabel('frequency[Hz]');
legend();

subplot(2,1,2)
hold on;
plot(f,imag(fftshift(H12)),"linewidth", 1.5,"LineStyle","--", "Color","r");
hold off;
xlabel('frequency[Hz]');
legend();



figure(5)
subplot(2,1,1)
hold on;
plot(f,abs(fftshift(H12)), "Color", "b");
hold off;
xlabel('frequency[Hz]');
legend();

subplot(2,1,2)
hold on;
% plot(f,180/pi*unwrap(angle(fftshift(H12))), ...
%     "Color", "b");
plot(f,180/pi*angle(fftshift(H12)), ...
    "Color", "b");
hold off;
xlabel('frequency[Hz]');
legend();




%% try with another function which is not a delta in freq domain




