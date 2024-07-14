function [H] = Fast_Fourier(IRxc,fs)

%% -----------------------------------------------------------------------
% INPUT
% IRxc = Struct with cross correlated impulse response data
%
% OUTPUT
% H = Matrix with transfer functions between mic 1-0, mic 2-0, mic 3-0, mic
% 4-0
%% -----------------------------------------------------------------------
% Data unwrapping
mic0 = IRxc.mic0; mic1 = IRxc.mic1;
mic2 = IRxc.mic2; mic3 = IRxc.mic3;
mic4 = IRxc.mic4;

% Determining shortest data set
% It's interesting here to use same vector length NFFT for every microphone
% data because the fft function needs an amplitude scaling that depends
% on the vevtor length NFFT. When they are the same size, this scalling is
% cancelled by the ratio of transfer function H.
NFFT = min([length(mic0), length(mic1), length(mic2), length(mic3),...
    length(mic4)]);


%Fast Fourier Transform
FFT0 = fft(mic0,NFFT);
FFT1 = fft(mic1,NFFT);
FFT2 = fft(mic2,NFFT);
FFT3 = fft(mic3,NFFT);
FFT4 = fft(mic4,NFFT);


% plot abs fft of microphones
% L=fs; %length of time signal
% fp=(0:L-1)/L*fs;
% figure(400+round(10*rand()));
% hold on;
% grid on;
% plot(fp(1:NFFT),abs(2/L*FFT0 ));
% plot(fp(1:NFFT),abs(2/L*FFT1 ));
% plot(fp(1:NFFT),abs(2/L*FFT2 ));
% plot(fp(1:NFFT),abs(2/L*FFT3 ));
% plot(fp(1:NFFT),abs(2/L*FFT4 ));
% xlim([0 5000]);
% xlabel('frequency[Hz]');
% ylabel('abs 2 fft/L')
% legend('0','1','2','3','4');
% title('fft of mics');
   
countnoise=0;
% compute transfer functions
for k=1:NFFT
   if(abs(FFT0(k))>1e-8) % if denominator isnt noise
    %Determination of transfer functions mic 0 - mic ii
    H(k,1) = FFT1(k)/FFT0(k); % transfer function of mic 1 compared to mic0
    H(k,2) = FFT2(k)/FFT0(k); % transfer fucntion of mic 2 compared to mic0
    H(k,3) = FFT3(k)/FFT0(k); % transfer function of mic 3 compared to mic0
    H(k,4) = FFT4(k)/FFT0(k); % transfer function of mic 4 compared to mic0
   else
    countnoise=countnoise+1;        
    H(k,1) = FFT1(k); % transfer function of mic 1 compared to mic0
    H(k,2) = FFT2(k); % transfer fucntion of mic 2 compared to mic0
    H(k,3) = FFT3(k); % transfer function of mic 3 compared to mic0
    H(k,4) = FFT4(k); % transfer function of mic 4 compared to mic0

    % % if you want the data to be nan
    % % H(k,1) = nan; % transfer function of mic 1 compared to mic0
    % % H(k,2) = nan; % transfer fucntion of mic 2 compared to mic0
    % % H(k,3) = nan; % transfer function of mic 3 compared to mic0
    % % H(k,4) = nan; % transfer function of mic 4 compared to mic0
    %     end
   end
   
end

fprintf(['Fast_Fourier called\n      number of noise detected: %d out' ...
         ' of %d NFFT data points\n\n'], countnoise, NFFT);

