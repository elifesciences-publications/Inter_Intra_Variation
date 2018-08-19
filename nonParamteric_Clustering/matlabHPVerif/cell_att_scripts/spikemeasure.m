function spikestats = spikemeasure(fs, vtrace, ctrace)
% Measure and return average values for spike threshold, spike width,
% spike height, spike half duration and rheobase. Requires voltage trace
% in mV and injected current ramp in pA.

% Last edited by Hugh Pastoll on 17-04-2013

% For code development
% fs = 20000;
% vtrace = c002_Membrane_Voltage_2*1e3;
% ctrace = c003_Current_2*1e12;


% Convert to double
type = whos('ctrace');
if strcmp(type.class,'cell')
    vtrace = cell2mat(vtrace);
    ctrace = cell2mat(ctrace);
end
vtrace = double(vtrace);
ctrace = double(ctrace);

% Indices of spike peaks
[spikemaxima, spikelocs] = findpeaks(vtrace, 'minpeakheight', 0);

%---------------------------------------------------
% Calculate spike threshold
%---------------------------------------------------

% (arbitrary) Threshold value
thresh_slope_val = 1;

% Only analyse up to the first 5 spikes - fast spiking at non-physiological
% rates can artificially raise the spike threshold.
max_spikes = 5;
% If there are fewer than 'max_spikes' in the trace, use the lower number
n_spike_limit = min([max_spikes length(spikelocs)]);

% Extract segment of trace prior to spike to find at which point the slope
% exceeds the (arbitrary) threshold value
segtime = 0.002; % 2 ms of trace prior to spike (in seconds)
segpts = segtime*fs;

% Store threshold values here
thresh_vm = zeros(length(spikelocs), 1);

% Optionally store extracted spikes trace segments and their slopes
% segs = zeros(length(spikelocs), segpts);
% slopes = zeros(length(spikelocs), segpts - 1);

% Calculate threshold for each spike
for i = 1:length(spikelocs)
    
    % Extract segment of trace prior to spike peak
    seg = vtrace(spikelocs(i) - segpts+1:spikelocs(i));
    
    % Optionally store segment
    %     segs(i, :) = vtrace(spikelocs(i) - segpts+1:spikelocs(i));
    
    % The slope between point in each segment
    segslopes = diff(seg);
    
    % Optionally store slopes
    %     slopes(i, :) = segslopes;
    
    % Remove all points in segment before last negative slope point to
    % ensure only points continuous with the spike are considered
    last_slope_decr_ind = find(segslopes<0, 1, 'last');
    segslopes(1:last_slope_decr_ind) = NaN;
    
    % Find the point in the segment where the slope first exceeds the
    % threshold slope value
    first_slope_thresh_ind = find(segslopes > thresh_slope_val, 1, 'first');
    if isempty(first_slope_thresh_ind)
        thresh_vm(i) = NaN;
    else
        thresh_vm(i) = seg(first_slope_thresh_ind);
    end
    
end

% Plot stored spikes and slopes

% figure(1)
% plot(segs')
% xlim([30 40])
% figure(2)
% plot(slopes')
% xlim([30 39])

% Store mean threshold value
spikestats(1) = nanmean(thresh_vm(1:n_spike_limit));

%---------------------------------------------------
% Calculate spike maximum
%---------------------------------------------------

% Store mean spike maximum
spikestats(2) = nanmean(spikemaxima(1:n_spike_limit));


%---------------------------------------------------
% Calculate spike halfwidth
%---------------------------------------------------

% Average half height
halfheight = ((spikestats(2)-spikestats(1))/2)+spikestats(1);

% Store width at half height values
spikehalfdur = zeros(n_spike_limit, 1);


% Calculate halfwidth for each spike
for i = 1:n_spike_limit
    
    % Extract segment of trace around spike peak
    seghalf = vtrace(spikelocs(i) - round(segpts/2+1):spikelocs(i)+round(segpts/2));
    
    % set all values less than halfheight to NaNs
    seghalf(seghalf > halfheight) = NaN;
    
    % get indexes of NaN values
    spikehighvals = isnan(seghalf);
    
    % find the index the first downward crossing of halfheight
    upindex = find(spikehighvals > 0, 1, 'first');
    
    % find the index of the last upward crossing of halfheight
    downindex = find(spikehighvals > 0, 1, 'last');
    
    spikehalfdur(i) = (downindex - upindex);
    
end

% Mean spike width at half height (in ms)
spikestats(3) = mean(spikehalfdur(1:n_spike_limit))/fs*1e3;

%---------------------------------------------------
% Calculate rheobase
%---------------------------------------------------

% Magnitude of injected current at time of first spike (in pA)
spikestats(4) = ctrace(spikelocs(1));

% figure(3)
% plotyy(c001_Time, c002_Membrane_Voltage_2, c001_Time, ctrace)


%---------------------------------------------------
% AHP measurements
%---------------------------------------------------

% AHP minimum (only consider first spike - assumes there is more than one
% spike)
if (length(spikelocs) < 2)
    warning('Too few spikes in ramp to calculate AHP')
end
spikestats(5) = min(vtrace(spikelocs(1):spikelocs(2)));

% FI slope. NB. Can produce NaNs if there are too few spikes in the trace
% to estimate the slope properly
% spikestats(6) = getFISlope(fs, vtrace, ctrace);

end








