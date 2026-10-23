function y = cwl_sosfiltfilt(sos, g, x)
% CWL_SOSFILTFILT  Zero-phase filtering of the columns of x with a second-order-section
% filter (sos, g). Since R2024b, filtfilt(sos, g, x) warns on every call that the inputs
% could also be read as a transfer function; they are read as SOS, which is what we
% want, so that warning is silenced here.
ws = warning('off', 'all');
restore = onCleanup(@() warning(ws));
y = filtfilt(sos, g, x);
end
