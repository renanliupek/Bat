function [IRxc] = Cross_correlation(IR, N, x, constants, fs)

%% -----------------------------------------------------------------------
% INPUT
% IR = Matrix with first second of impuls response data
% N = Amount of datapoints for first sound to travel through the tube
% x = Distance between mic 0 and mic 1-4
% constants = Struct with all used constants (rho, c, T, P)
% fs = Amount of datapoints per second
%
% OUTPUT
% IRxc = struct with cross correlated impuls response of all mics
%% -----------------------------------------------------------------------
% Cross correlation of audio data between ref mic (mic 0) and other
% mics(1-4)
for ii = 2:length(IR(1,:))
    C = xcorr(IR(1:N+1, 1), IR(1:N+1, ii));
    [~, I(ii-1)] = max(C, [], 1);
end

% Difference between the first peak in audio data between ref mic and other
% mics
diff_real = N*ones(size(I)) - I;

% Expected difference between first peak in audio data between ref and
% other mics
diff_expected = round((x./constants.c).*fs);

% Number of samples the impulse response should be shifted forward to
% correlate with mic 0
s = diff_expected - diff_real;

shift = [s(2)+s(3)+s(4), s(1)+s(3)+s(4), s(1)+s(2)+s(4), s(1)+s(2)+s(3)];

q = min(shift);

data_cut = shift-q;

IRxc.mic0 = IR(:,1);
IRxc.mic1 = IR(data_cut(1)+1:end, 2);
IRxc.mic2 = IR(data_cut(2)+1:end, 3);
IRxc.mic3 = IR(data_cut(3)+1:end, 4);
IRxc.mic4 = IR(data_cut(4)+1:end, 5);


end
