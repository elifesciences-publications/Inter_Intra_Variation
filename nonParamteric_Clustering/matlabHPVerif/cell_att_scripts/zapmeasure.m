function zap_vals = zapmeasure(fs, vtrace, ctrace)
% Computes the resonant frequency (f_max) and resonance magnitude (q_val)
% by deriving an impedance amplitude profile (ZAP) from voltage trace and
% injected current with linear chirp waveform ranging from 1 - 20 Hz.
% Requires voltage trace in mV and current injection waveform in pA.

% Last edited by Hugh Pastoll 17-04-2013

% % For code development
% vtrace = c002_Membrane_Voltage_2*1e3;
% ctrace = c003_Current_2*1e12;
% fs = 20000;
% Convert to Double
% Convert to double
type = whos('ctrace');
if strcmp(type.class,'cell')
    vtrace = cell2mat(vtrace);
    ctrace = cell2mat(ctrace);
end
vtrace = double(vtrace);
ctrace = double(ctrace);
% FFT points and frequencies
fftpts = length(vtrace);
freqs = fs*(0:fftpts/2-1)/fftpts; % frequency range (for plotting)

% FFTs of voltage and current traces
y_v = fft(vtrace, fftpts);
y_c = fft(ctrace, fftpts);

% Amplitude spectra
pyy_v = abs(y_v)/(fftpts/2);
pyy_c = abs(y_c)/(fftpts/2);


% First half of 'half' of frequency spectra
yrange_v = pyy_v(1:fftpts/2);
yrange_c = pyy_c(1:fftpts/2);

% Ratio of spectra
fft_ratio = yrange_v./yrange_c;

% Smooth the fft ratio
smoothing = 7;
fft_ratio_smooth = filter(1/smoothing*ones(smoothing,1),1,fft_ratio);

% Minima and maxima of frequency range
freq_min = 1;
freq_max = 20;

% Indices for extracting segments of traces between frequency minima and
% maxima
[~, freq_min_ind] = min(abs(freqs-freq_min));
[~, freq_max_ind] = min(abs(freqs-freq_max));

% Extract segments to calculate F_max and Q
freqs_seg = freqs(freq_min_ind:freq_max_ind);
fft_ratio_seg = fft_ratio_smooth(freq_min_ind:freq_max_ind);

% F_max and Q values
[~, f_max_ind] = max(fft_ratio_seg);
zap_vals(1) = freqs_seg(f_max_ind);
zap_vals(2) = fft_ratio_seg(f_max_ind)./fft_ratio_seg(1);


% figure(1)
% m = 3;n = 1;
% 
% % Plot amplitude spectrum of voltage
% subplot(m, n, 1)
% semilogx(freqs, yrange_v)
% xlim([freq_min freq_max])
% 
% % Plot amplitude spectrum of injected current
% subplot(m, n, 2)
% semilogx(freqs, yrange_c)
% xlim([freq_min freq_max])
% 
% % Plot amplitude spectra ratio
% subplot(m, n, 3)
% plot(freqs, fft_ratio_smooth)
% xlim([freq_min freq_max])

end