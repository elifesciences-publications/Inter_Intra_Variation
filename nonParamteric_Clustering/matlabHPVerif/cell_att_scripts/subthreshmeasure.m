function subthreshstats = subthreshmeasure(fs, vtraces_raw, ctraces)
warning('off','MATLAB:colon:nonIntegerIndex');
opts = optimoptions('lsqcurvefit','Display','off');
% Do not analyse the trace with no current injection
zero_step_trace_ind = 3;
% Convert to double
type = whos('ctraces');
if strcmp(type.class,'cell')
    vtraces_raw = cell2mat(vtraces_raw);
    ctraces = cell2mat(ctraces);
end
% Smooth traces
smoothing = 5;
vtraces = filter(1/smoothing*ones(smoothing,1),1,vtraces_raw);

%---------------------------------------------------------------
% Calculate resting membrane potential
%---------------------------------------------------------------

% Mean of the first second of all traces to get resting potential
rest_end = 1; % in seconds

% near but not at beginning of trace to avoid smoothing artefact
rest_start_ind = 5;
% end of 1st second of trace
rest_end_ind = rest_end*fs;

% Return mean resting membrane potential
subthreshstats(1) = mean(mean(vtraces(rest_start_ind:rest_end_ind, :)));


%---------------------------------------------------------------
% Calculate input resistance
%---------------------------------------------------------------


% Use final two seconds during current steps to calculate steady state
ss_start = 2; % in seconds
ss_dur = 2; % in seconds

ss_start_ind = ss_start*fs+1; % beginning of steady state
ss_end_ind = (ss_start+ss_dur)*fs; % end of steady state

% Amplitude of the voltage deflection
delta_v = mean(vtraces(ss_start_ind:ss_end_ind, :))...
    - mean(vtraces(rest_start_ind:rest_end_ind, :));

% Amplitude of the current step
delta_c = mean(ctraces(ss_start_ind:ss_end_ind, :))...
    - mean(ctraces(rest_start_ind:rest_end_ind, :));

% Input resistance R = V/I
input_resistance = delta_v./delta_c*1e3; % in M? (G? without conversion)

% Remove zero step trace
input_resistance(zero_step_trace_ind) = [];

% Return mean input resistance value of negative steps
subthreshstats(2) = mean(input_resistance(1:2));

%---------------------------------------------------------------
% Calculate sag coefficient
%---------------------------------------------------------------


% Extract 75 ms segment of trace 2 ms after current step onset to find sag
% maxima and minima

sag_start = rest_end + 0.002; % in seconds
sag_extract_dur = 0.075; % in seconds

sag_start_ind = round(sag_start*fs);
sag_end_ind = round((sag_start+sag_extract_dur)*fs);

% Extract portions containing sag and normalise to to resting potential.
% Then take absolute value to permit taking the maximum value of all traces
sag_extract = abs(bsxfun(@minus, vtraces(sag_start_ind:sag_end_ind, :),...
    mean(vtraces(rest_start_ind:rest_end_ind, :))));

% Sag is ratio of absolute value of the steady state voltage deflection to
% maximum voltage deflection
sag_vals = abs(delta_v)./max(sag_extract);

% Remove zero step trace
sag_vals(zero_step_trace_ind) = [];

% Return mean sag coefficient
subthreshstats(3) = mean(sag_vals);


%---------------------------------------------------------------
% Calculate membrane time constant
%---------------------------------------------------------------


% Extract 12 ms segment of trace 2 ms after current step onset to find time
% constant

tau_start = rest_end + 0.002; % in seconds
tau_dur = 0.012; % in seconds

tau_start_ind = round(tau_start*fs);
tau_end_ind = round((tau_start+tau_dur)*fs);

% Extract portions containing sag and normalise to to resting potential.
% Then take absolute value to permit taking the maximum value of all traces
tau_extract = abs(bsxfun(@minus, vtraces(tau_start_ind:tau_end_ind, :),...
    mean(vtraces(rest_start_ind:rest_end_ind, :))));

% Remove zero step trace
tau_extract(:, zero_step_trace_ind) = [];


% Get timebase of voltage deflections (in ms)
tau_x = (1:size(tau_extract, 1))/fs*1000;

% Number of traces to fit
ntraces = size(tau_extract, 2);

% Store time constant fits
tau = zeros(ntraces, 1);

% Store fits
tau_hat = zeros(size(tau_extract));


% Fit time constant to all traces
for n = 1:ntraces
    
    % Use anonymous function with voltage offset, amplitude and time constant
    % parameters to estimate tau. 'h' is a function handle
    h = @(m, x) m(1) + m(2) * exp(x./m(3));
    
    % Call fit function with anonymous function handle, initial guesses at
    % parameters, the timebase of voltage deflections and the trace.
    
    [parameters] = lsqcurvefit(h, [0 1 -10].', tau_x.',double(tau_extract(:, n)),[0 -20 -50],[20 0 0],opts);
    
    % Membrane time constant (in ms)
    tau(n) = -parameters(3);
    
    
    % Calculate fit
    tau_hat(:, n) = parameters(1) + ...
        parameters(2)*exp(tau_x./parameters(3));
    
end


 
% % Plot in the same set of axes
% plot(tau_x, tau_extract, tau_x, tau_hat);

% % Optionally allow user to check that fits are accurate
% pause;

% Return median membrane time constant (averages two 'middle' values)
subthreshstats(4) = median(tau);

end



