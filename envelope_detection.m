function [signal_t, signal_f] = envelope_detection(signal, fs, fs_2)

envelope_tc = abs(hilbert(signal));
signal_t = resample(envelope_tc, fs, fs_2);
signal_f = fftshift(fft(signal_t));

end