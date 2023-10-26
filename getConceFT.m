

% Add path to the ConceFT codes in your environment
addpath('/Volumes/Riki2tb/PRUV/Sleep/Time-Frequency-Analysis-Matlab-code');
addpath('/Volumes/Riki2tb/PRUV/Sleep/Time-Frequency-Analysis-Matlab-code/tool') ;
addpath('/Volumes/Riki2tb/PRUV/Sleep/Time-Frequency-Analysis-Matlab-code/Morse') ;



base_dir = "/Volumes/Riki2tb/PRUV/Sleep";
data_dir = fullfile(base_dir, "Dream_Spindle_Data");

subject = 1;
subject = int2str(subject);
filename = "excerpt" + subject + ".edf";
filename_hyp = "Hypnogram_excerpt" + subject + ".txt";
filename_expert = "Visual_scoring1_excerpt" + subject + ".txt";


% the sampling rate for the simulated signal
Hz = 50;





%%% Get the excerpt data

if subject == "3"

    filepath = fullfile(data_dir, filename);
    data = edfread(filepath);

    cz_temp = data(:, "C3_A1");
    cz_temp = timetable2table(cz_temp);
    cz_temp = table2array(cz_temp(:, 2));

%%% This is just for subject 3
    cz_temp2 = zeros(height(cz_temp)*100, 1);
    for i = 0:height(cz_temp)-1
        cz_temp2(i*100+1:(i+1)*100, :) = cell2mat(cz_temp(i+1, :));
    end
    
    % Upsample the signal
    cz = cz_temp2;
    

elseif subject == "1"
    
    filename = "c301.mat";
    filepath = fullfile(data_dir, filename);
    cz = (load(filepath).c3)';
    cz = decimate(cz, 2);

else
    filepath = fullfile(data_dir, filename);
    data = edfread(filepath);

    cz_temp = data(:, "CZ_A1");
    cz_temp = timetable2table(cz_temp);
    cz_temp = table2array(cz_temp(:, 2));

    cz_temp2 = zeros(height(cz_temp)*400, 1);
    for i = 0:height(cz_temp)-1
        cz_temp2(i*400+1:(i+1)*400, :) = cell2mat(cz_temp(i+1, :));
    end
    cz = decimate(cz_temp2, 4);
end



J = 2;
Q = 20;
FrequencyAxisResolution = 1e-2;


[tfrsqtic, ConceFT] = getConceFT2(cz, J, Q, FrequencyAxisResolution, Hz);


function [tfrsqtic, ConceFT] = getConceFT2(signal, J, Q, FrequencyAxisResolution, Hz)
%% setup parameters for the SST or ConceFT
% signal: 1D array [# of time points, 1]
% This is the signal that we want to calculate ConceFT repreentation for.
% J: int
% number of chosen orthonormal windows for ConceFT
% Q: int
% number of random linear combinations of chosen windows
% FrequencyAxisResolution: float 
% The resolution of frequency axis. 1e-3 or 1e-4 is recommended.
% Hz: float
% The sampling frequency of signals



%% Other parameters (should not change)

% the window length. Ideally, it should be chosen so that
% roughly 7-10 oscillations (ignore the multiples) are
% included in the window.
factor = 2;
WindowLength = factor * Hz+1 ; % 
% this is the bandwith of the chosen window. See hermf.m
% in the attached code for details.
WindowBandwidth = 10;
SamplingRate = Hz;


% Setup the frequency range for the analysis
% The allowed range is 0-0.5
% This part might be tricky. This 0-0.5 limitation is
% setup under the assumption that the sampling rate is 1Hz
% After the analysis, the sampling rate should be adjusted
% so that the result is associated with the original
% sampling rate.
% In this example, the true range is [0, 0.5]*SamplingRate
HighFrequencyLimit = 20/Hz; % Get frequency axis upto 20 Hz
LowFrequencyLimit = 0;


% the frequency axis resolution in the final time-frequency representation
HOP = 1;

% call the main code, which is the ConceFT based on
% synchrosqueezed short time Fourier transform (STFT)
% Output:
% tfr: STFT result
% tfrtic: frequency axis tic for the STFT
% tfrsq: synchrosqueezed STFT (it is equivalent to running ConceFT only one time)
% ConceFT: ConceFT of synchrosqueezed STFT.
% tfrsqtic: frequency axis tic for the tfrsq and ConceFT
[tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_sqSTFT_C(signal, LowFrequencyLimit, ...
    HighFrequencyLimit, FrequencyAxisResolution, HOP, WindowLength, ...
    J, WindowBandwidth, Q, 1, 0);
end

