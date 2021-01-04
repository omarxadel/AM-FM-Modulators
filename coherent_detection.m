function [signal_t, signal_f] = coherent_detection(signal, carrier, fs, fs_2, bw)

B = signal.*carrier;
Fnormalized = fs/2;
cutoffFreq = bw/Fnormalized;
[zeroes, poles] = butter(5, cutoffFreq, 'low');
B = filtfilt(zeroes, poles, B);
signal_t = resample(B, fs, fs_2);
signal_f = fftshift(fft(signal_t));

end