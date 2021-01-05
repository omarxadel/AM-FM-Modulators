function [signal_t, signal_f] = envelope_detection(signal, fs, fs_2)

envelope = abs(hilbert(signal));
signal_t = resample(envelope, fs, fs_2);
signal_f = fftshift(fft(signal_t));

end