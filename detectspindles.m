
function [precision, sensitivity, f1, idx_detect, dur_stats, density_stats] = detectspindles(ConceFT, tfrsqtic, epsilon, delta, fmin, fmax, tmin, tmax, min_distance_ms, FrequencyAxisResolution, idx_expert, hypnogram, Hz, range)

%%% Parameters

% ConceFT: 2D array
%   ConceFT representation of the signal. 2d matrix [# of freq index, # of
%   time points]. Use the output of getConceFT.m
% tfrsqtic: 1D array
%   frequenxy axis of the ConceFT normalized by sampling frequency. Use the
%   output of getConceFT.m
% epsilon: float
%   a threshold parameter defined in the paper, (0, 1)
% delta: float
%   a threshold parameter defined in the paper, (0, Infinity)
% fmin, fmax: float
%   used to specify the sigma frequency band range. In the paper, fmin=12, fmax=15
% tmin, tmax: float
%   used to specify the minimum and maximum length of spindles (cutoffs in milliseconds). In
%   the paper, tmin = 300, tmax = 3000
% min_distance_ms: float
%   if the gap between neighboring spindles are less than min_distance_ms,
%   they will be merged. In the paper, min_ditance_ms = 3000
% FrequencyAxisResolution: float
%   FrequencyAxisResolution is the same parameter used to obtain ConceFT. 
% idx_expert: 1d array, column vector [# of time points, 1]
%   array that indicates whether a spindle exist at each time point, labeled
%   by expert. Should be the same length as # of time points of ConceFT
% hypnogram: 1d array, column vector [# of time points, 1]
%   array that indicates the sleep stage at each time point, labeled.
%   by expert. Should be the same length as # of time points of ConceFT.
%   For the N2 sleep stage, hypnoram entry must be 2
% Hz: float
%   frequency of the signals
% range: float
%   frequency range used to calculate sigma band energy. The code averages 
%   sigma band power around the extracted curve +/- range
%   
%
%%%



%% Extract Curves from the specified frequency range
ConceFT_ftic = Hz * tfrsqtic;
freq_index = find(ConceFT_ftic >= fmin & ConceFT_ftic <= fmax);
valid_freq = ConceFT_ftic(freq_index);
ConceFT_window = ConceFT(freq_index, :);
[c_CF] = CurveExt_M(abs(ConceFT_window)', .001);
c_CF_freq = valid_freq(c_CF);

%Calculate the energy band around the extracted curve (ConceFT)
curve_index = c_CF + max(find(ConceFT_ftic < fmin));
energy = zeros(length(curve_index), 1);
Q = round(range/(Hz*FrequencyAxisResolution)); % express the range in terms of number of indices

for jj = 1: length(curve_index)
    energy(jj) = sum(ConceFT(curve_index(jj)-Q:curve_index(jj)+Q, jj).^2)/(range*2); % This is the sigma band energy density
end


%% Calculate the delta/theta/alpha band power density
idx_delta = find(ConceFT_ftic >0.5 & ConceFT_ftic <= 4);
idx_theta = find(ConceFT_ftic >4 & ConceFT_ftic <= 8);
idx_alpha = find(ConceFT_ftic >8 & ConceFT_ftic <= fmin);

delta_amp = ConceFT(idx_delta, :);
delta_power = delta_amp.^2;
delta_power = sum(delta_power, 1)/(4-0.5);

theta_amp = ConceFT(idx_theta, :);
theta_power = theta_amp.^2;
theta_power = sum(theta_power, 1)/(8-4);

alpha_amp = ConceFT(idx_alpha, :);
alpha_power = alpha_amp.^2;
alpha_power = sum(alpha_power, 1)/(fmin-8);

%% Apply the threshold to normalized sigma power
normalized_sigma_power = energy' ./(delta_power + theta_power + alpha_power + energy');
idx_sigma_normal = find(normalized_sigma_power >= epsilon);



%% Get the indeices for thresholded amlitude in sigma band
idx_sigma2 = find(ConceFT_ftic >fmin & ConceFT_ftic <= fmax);
sigma_amp2 = ConceFT(idx_sigma2, :);
amplitude = sum(sigma_amp2, 1);

%% Define hard and soft thresholds
hard_thr = mean(amplitude) + delta * std(amplitude);
soft_thr = 0.5 * hard_thr;

idx_hard = find(amplitude > hard_thr);
idx_soft = find(amplitude > soft_thr);
idx_zc_soft = events_to_index_flatten(idx_soft);

spindles = process_spindles(idx_hard, idx_sigma_normal, idx_zc_soft, min_distance_ms, Hz, tmin, tmax);


%% Create the 1D array, indicating whether at each time point it was predicted to be a spindle
idx_detect = zeros(length(idx_expert), 1);
if numel(spindles) ~= 0
    for i = 1:size(spindles, 1)
        idx_detect(spindles(i, 1):spindles(i, 2)) = 1;
    end
end


%% Now calculate Precision, Sensitivity, F1 score of N2 stage
% Only interested in N2 sleep stage
idx_expert = idx_expert(hypnogram == 2);
idx_detect_n2 = idx_detect(hypnogram == 2);

% Get the spindle detection accuracy by event
idx_expert_inp = find(idx_expert == 1)';
idx_detect_inp = find(idx_detect_n2 == 1)';
spindles_detected = events_to_index(idx_detect_inp);
spindles_expert = events_to_index(idx_expert_inp);

overlap_thresholds = [0.2];
[precision, sensitivity, f1, overlap_valid] = f1_scores(spindles_detected, spindles_expert, overlap_thresholds);


%% Calculate durations of spindles by conceft-s/expert on N2

tp_gs_index = any(overlap_valid, 1);
tp_detected_index = any(overlap_valid, 2);

tp_gs = spindles_expert(tp_gs_index, :);
tp_detected = spindles_detected(tp_detected_index, :);

durs_detected = (tp_detected(:, 2) - tp_detected(:, 1))/Hz;
durs_expert = (tp_gs(:, 2) - tp_gs(:, 1))/Hz;

dur_stats = cell(2, 1);
dur_stats{1} = durs_detected;
dur_stats{2} = durs_expert;

%% This section for spindle density calculation for n2 expert, n2 conceft
n2_expert_density = size(spindles_expert, 1)/(length(find(hypnogram == 2))/Hz);
n2_conceft_density = size(spindles_detected, 1)/(length(find(hypnogram == 2))/Hz);
density_stats = [n2_expert_density, n2_conceft_density];

end





function overlap = get_overlap(spindles_detected, spindles_gs)
n_detected_spindles = size(spindles_detected, 1);
n_gs_spindles = size(spindles_gs, 1);
% The (relative) overlap between each pair of detected spindle and gs spindle
overlap = zeros(n_detected_spindles, n_gs_spindles);
for i = 1:n_detected_spindles
    for j = 1:n_gs_spindles
        idx_detected = spindles_detected(i, :);
        idx_gs = spindles_gs(j, :);

        % [start, stop) indices of the detected spindle and of the gs spindle
        idx_range_detected = idx_detected(1):idx_detected(2);
        idx_range_gs = idx_gs(1):idx_gs(2);

        % Calculate intersect and union of the spindle indices
        intersect = intersect1d(idx_range_detected, idx_range_gs);
        union = union1d(idx_range_detected, idx_range_gs);

        % Overlap of a detected spindle and a gs spindle is defined as the intersect over the union
        overlap(i, j) = length(intersect) / length(union);
    end
end
end

function result = intersect1d(a, b)
[~, ia, ~] = intersect(a, b);
result = a(ia);
end

function result = union1d(a, b)
result = unique([a, b]);
end



function [n_true_positives, overlap_valid] = get_true_positives(spindles_detected, spindles_gs, overlap_thresholds)
% If either there is no spindle detected or the gold standard doesn't contain any spindles, there can't be any true
% positives
if isempty(spindles_detected) || isempty(spindles_gs)
    n_true_positives = zeros(size(overlap_thresholds));
    return;
end

% Get the overlaps in format (n_detected_spindles, n_gs_spindles)
overlap = get_overlap(spindles_detected, spindles_gs);
% Make sure there is at max one detection per gs event
[~, maxIndicesCols] = max(overlap, [], 1);
linearIndicesCols = sub2ind(size(overlap), maxIndicesCols, 1:size(overlap,2));
overlap_valid1 = zeros(size(overlap));
overlap_valid1(linearIndicesCols) = overlap(linearIndicesCols);

% Make sure there is at max one gs event per detection
[~, maxIndicesRows] = max(overlap_valid1, [], 2);
linearIndicesRows = sub2ind(size(overlap_valid1), 1:size(overlap_valid1,1), maxIndicesRows');
overlap_valid = zeros(size(overlap_valid1));
overlap_valid(linearIndicesRows) = overlap_valid1(linearIndicesRows);

n_true_positives = zeros(size(overlap_thresholds));

% Calculate the valid matches (true positives) depending on the overlap threshold
for idx = 1:numel(overlap_thresholds)
    % All remaining values > overlap_threshold are valid matches (true positives)
    matches = find(overlap_valid > overlap_thresholds(idx));
    n_true_positives(idx) = numel(matches);
end
end



function [precision, recall, f1, overlap_valid] = f1_scores(spindles_detected, spindles_gs, overlap_thresholds)
n_detected_spindles = size(spindles_detected, 1);
n_gs_spindles = size(spindles_gs, 1);

% Get the number of true positives per overlap threshold
[n_true_positives, overlap_valid] = get_true_positives(spindles_detected, spindles_gs, overlap_thresholds);

[precision, recall, f1] = metric_scores(n_detected_spindles, n_gs_spindles, n_true_positives);
end


function [precision, recall, f1] = metric_scores(n_detected_spindles, n_gs_spindles, n_true_positives)
if (n_detected_spindles == 0) && (n_gs_spindles == 0)
    % If there are no spindles detected and the gold standard doesn't contain any spindles, the precision, recall,
    % and f1 score are defined as one
    precision = ones(size(n_true_positives));
    recall = ones(size(n_true_positives));
    f1 = ones(size(n_true_positives));
elseif (n_detected_spindles == 0) || (n_gs_spindles == 0)
    % If either there are no spindles detected or the gold standard doesn't contain any spindles, there can't be any
    % true positives, and precision, recall, and f1 score are defined as zero
    precision = zeros(size(n_true_positives));
    recall = zeros(size(n_true_positives));
    f1 = zeros(size(n_true_positives));
else
    % Precision is defined as TP/(TP+FP)
    precision = n_true_positives ./ n_detected_spindles;
    % Recall is defined as TP/(TP+FN)
    recall = n_true_positives ./ n_gs_spindles;
    % f1 score is defined as the harmonic mean between precision and recall
    f1 = harmmean([precision, recall], 2);
end
end



function xidx = events_to_index(x)
% Convert a continuous vector of indices into a 2D array (start, end).
% Split indices where it stopped :
sp = splitVector(x, find(diff(x) ~= 1));
% Return (start, end) :
xidx_temp = cellfun(@(k) [k(1), k(end)], sp, 'UniformOutput', false);
xidx = zeros(numel(xidx_temp), 2);
for i = 1:numel(xidx_temp)
    xidx(i, :) = cell2mat(xidx_temp(i));
end
%xidx = cell2mat(xidx);
end

function xidx = events_to_index_flatten(x)
% Convert a continuous vector of indices into a 1D array (start, end).
% Split indices where it stopped :
sp = splitVector(x, find(diff(x) ~= 1));
% Return (start, end) :
xidx_temp = cellfun(@(k) [k(1), k(end)], sp, 'UniformOutput', false);
xidx = cell2mat(xidx_temp);
end



function result = splitVector(vec, indices)
% Split a vector into cells at specified indices.
if any(indices==1) && any(indices==numel(vec))
elseif any(indices==1)
    indices = [indices, numel(vec)];
elseif any(indices==numel(vec))
    indices = [1, indices];
else
    indices = [1, indices, numel(vec)];
end
result = arrayfun(@(i) vec(indices(i)+1:indices(i+1)), 1:numel(indices)-1, 'UniformOutput', false);
end


function f_index = events_distance_fill(index, min_distance_ms, sf)
% Remove events that do not have the good duration.

% Convert min_distance_ms
min_distance = min_distance_ms / 1000 * sf;
idx_diff = diff(index);
condition = idx_diff > 1;
idx_distance = find(condition);
distance = idx_diff(condition);
bad = idx_distance(distance < min_distance);

% Fill gap between events separated with less than min_distance_ms
if ~isempty(bad)
    fill = [];
    for j = 1:numel(bad)
        fill = [fill, index(bad(j)) + 1:index(bad(j) + 1) - 1];
    end
    f_index = sort([index, fill]);
else
    f_index = index;
end
end



function output = process_spindles(idx_hard, idx_sigma_normal, idx_zc_soft, min_distance_ms, sf, tmin, tmax)
% Process spindles based on the provided indices and parameters.
output = [];
if numel(idx_hard) > 0
    % Initialize spindles vector
    idx_spindles = [];

    % Keep only period with high relative sigma power
    [idx_hard, ~, ~] = intersect(idx_hard, idx_sigma_normal);


    % Fill gap between events separated by less than min_distance_ms
    idx_hard = events_distance_fill(idx_hard, min_distance_ms, sf);

    % Get where spindles start / end
    if isempty(idx_hard) || length(idx_hard) <= 1
        output = [];
        return;
    end

    idx = events_to_index(idx_hard);
    idx_start = idx(:, 1);
    idx_stop = idx(:, 2);

    % Find true beginning / end using soft threshold
    for s = idx_start'
        d = s - idx_zc_soft;
        % Find distance to nearest soft threshold crossing before start
        soft_beg = min(d(d > 0));
        % Find distance to nearest soft threshold crossing after end
        soft_end = min(abs(d(d < 0)));
        idx_spindles = [idx_spindles, (s - soft_beg):(s + soft_end)];
    end

    % Fill gap between events separated by less than min_distance_ms
    idx_spindles = events_distance_fill(idx_spindles, min_distance_ms, sf);

    if isempty(idx_spindles) || length(idx_spindles) <= 1
        output = [];
        return;
    end

    % Get duration
    idx = events_to_index(idx_spindles);
    idx_start = idx(:, 1);
    idx_stop = idx(:, 2);
    duration_ms = (idx_stop - idx_start) * (1000 / sf);

    % Remove events with bad duration
    good_dur = find(duration_ms > tmin & duration_ms < tmax);

    if isempty(idx_spindles)
        output = [];
        return;
    end

    output = [idx_start(good_dur), idx_stop(good_dur)];
    output = unique(output, 'rows');
end
end




