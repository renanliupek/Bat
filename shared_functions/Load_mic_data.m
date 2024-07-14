function [fs, IRmic] = Load_mic_data(Sample, Backing)

%% ------------------------------------------------------------------------
% INPUT
% Sample = File name of the sample
% Backing = Backing used in the measurement (Anechoic, Hard or Open)

% OUTPUT
% fs = Amount of measurement point per second
% IRmic = Matrix with full audio data of mic 0 (column 1), 1, 2, 3 and 4(column
% 5)
%% -----------------------------------------------------------------------
    % If-loop to read to files for the right backing
    if strcmp(Backing, 'Anechoic')
        File = ["IRmica0.wav", "IRmica1.wav", "IRmica2.wav", ...
            "IRmica3.wav", "IRmica4.wav"];
    elseif strcmp(Backing, 'Hard')
        File = ["IRmich0.wav", "IRmich1.wav", "IRmich2.wav", ...
            "IRmich3.wav", "IRmich4.wav"];
    elseif strcmp(Backing, 'Open')
        File = ["IRmico0.wav", "IRmico1.wav", "IRmico2.wav", ...
            "IRmico3.wav", "IRmico4.wav"];
    end
    
    % File loading for microphone 0, 1, 2, 3 & 4
    
    [IRmic(:,1),fs] = audioread(fullfile(Sample, File(1)));     
    [IRmic(:,2),~]  = audioread(fullfile(Sample, File(2)));
    [IRmic(:,3),~]  = audioread(fullfile(Sample, File(3)));
    [IRmic(:,4),~]  = audioread(fullfile(Sample, File(4)));    
    [IRmic(:,5),~]  = audioread(fullfile(Sample, File(5)));
%     disp(fs);
end
