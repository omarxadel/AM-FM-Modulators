% VARIABLES
fig_i = 1;

% 2.Frequency Modulation
% 2.1.Reading the audio input.
[signal,Fs] = audioread('eric.wav'); 

FrequencyDomainSignal = fftshift(fft(signal));
f = linspace(-Fs/2,Fs/2,length(FrequencyDomainSignal));
figure(fig_i);
fig_i = fig_i + 1;
plot(f, abs(FrequencyDomainSignal));
title('Spectrum of the input signal');
xlabel('Frequency'); 
ylabel('Value');

% 2.2.Filtering frequencies > 4KHz.
BW = 4000;
filter = ones(length(FrequencyDomainSignal), 1);
   
    for i = 1:length(FrequencyDomainSignal)
        if f(i)<-BW || f(i)>BW
            filter(i) = 0;
        end
    end
    
FilteredSignal_f = filter.*FrequencyDomainSignal;

figure(fig_i);
fig_i = fig_i + 1;
f_filtered = linspace(-Fs/2,Fs/2,length(real(FilteredSignal_f)));
plot(f_filtered, abs(FilteredSignal_f));
title('Filtered signal in frequency domain');
xlabel('Frequency'); 
ylabel('Value');

% 2.3.Filtered signal in time domain.
FilteredSignal_t = real(ifft(ifftshift(FilteredSignal_f)));

figure(fig_i);
fig_i = fig_i + 1;
t = linspace(0,length(FilteredSignal_t)/Fs, length(FilteredSignal_t));
plot(t, real(FilteredSignal_t));
title('Filtered signal in time domain');
xlabel('Time'); 
ylabel('Value');

% 2.4.Playing the filtered audio signal (sound).
sound(FilteredSignal_t, Fs);
pause(8);

% 2.5.Generation of the NBFM signal.
Fc = 100000;
resampleFrequency = 5*Fc;
resampledSignal = resample(FilteredSignal_t, resampleFrequency,Fs);

t = linspace(0,length(resampledSignal)/resampleFrequency, length(resampledSignal));
A = max(abs(resampledSignal));
kf = pi;
integrated_message = kf*cumsum(resampledSignal).';

FM_t = A.*cos((2*pi*Fc*t) + integrated_message);
FM_f = fftshift(fft(FM_t));

f = linspace(-Fs/2, Fs/2, length(FM_f));
figure(fig_i);
fig_i = fig_i + 1;
plot(f, FM_f);
title('Modulated signal');
xlabel('Frequency'); 
ylabel('Value'); 

% 2.6.Demodulation of the NBFM signal.
diff_FM = diff(FM_t);

FM_demod_t = resample(abs(hilbert(diff_FM)), Fs, resampleFrequency);

t = linspace(0, length(FM_demod_t)/Fs, length(FM_demod_t));
figure(fig_i);
fig_i = fig_i + 1;
plot(t, FM_demod_t)
title('Demodulated signal at time domain');
xlabel('Time'); 
ylabel('Value'); 
sound(FM_demod_t, Fs);
