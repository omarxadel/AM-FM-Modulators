% VARIABLES
fig_i = 1;
 
% 2.Frequency Modulation
% 2.1.Reading the audio input.
[signal,Fs] = audioread('eric.wav'); 
 
FrequencyDomainSignal = fftshift(fft(signal)); %Get the signal in frequency domain
f = linspace(-Fs/2,Fs/2,length(FrequencyDomainSignal)); %Frequency range to plot the spectrum of the signal
figure(fig_i);
fig_i = fig_i + 1;
plot(f, abs(FrequencyDomainSignal));
title('Spectrum of the input signal');
xlabel('Frequency'); 
ylabel('Value');
 
% 2.2.Filtering frequencies > 4KHz.
BW = 4000; %To obtain the desired band limited signal
%Ideal filter to remove frequencies greater than 4KHz
filter = ones(length(FrequencyDomainSignal), 1);
   
    for i = 1:length(FrequencyDomainSignal)
        if f(i)<-BW || f(i)>BW
            filter(i) = 0;
        end
    end
    
FilteredSignal_f = filter.*FrequencyDomainSignal;%Filtered signal in frequency domain
 
figure(fig_i);
fig_i = fig_i + 1;
f_filtered = linspace(-Fs/2,Fs/2,length(real(FilteredSignal_f))); %Frequency range to plot the filtered signal in frequency domain
plot(f_filtered, abs(FilteredSignal_f));
title('Filtered signal in frequency domain');
xlabel('Frequency'); 
ylabel('Value');
 
% 2.3.Filtered signal in time domain.
FilteredSignal_t = real(ifft(ifftshift(FilteredSignal_f))); %Inverse fourier to get the filtered signal in time domain
 
figure(fig_i);
fig_i = fig_i + 1;
t = linspace(0,length(FilteredSignal_t)/Fs, length(FilteredSignal_t)); %Time range to plot the the filtered signal in time domain
plot(t, real(FilteredSignal_t));
title('Filtered signal in time domain');
xlabel('Time'); 
ylabel('Value');
% 2.4.Playing the filtered audio signal (sound).
sound(FilteredSignal_t, Fs);
pause(8);
 
% 2.5.Generation of the NBFM signal.
Fc = 100000; %Carrier frequency
resampleFrequency = 5*Fc; %Sampling frequency
resampledSignal = resample(FilteredSignal_t, resampleFrequency,Fs); %Resample the filtered signal to the new sampling frequency instead of Fs which is the original sampling frequency of the signal (=48KHz)
 
t = linspace(0,length(resampledSignal)/resampleFrequency, length(resampledSignal)); %Time range of the new resampled signal
A = max(abs(resampledSignal)); %Peak of the resampled signal
kf = pi; %Frequency deviation constant
integrated_message = kf*cumsum(resampledSignal).'; 
 
FM_t = A.*cos((2*pi*Fc*t) + integrated_message); %Modulated signal in time domain
FM_f = fftshift(fft(FM_t)); %Modulated signal in frequency domain
 
f = linspace(-Fs/2, Fs/2, length(FM_f)); %Frequency range of the modulated signal in frequency domain
figure(fig_i);
fig_i = fig_i + 1;
plot(f, FM_f);
title('Modulated signal');
xlabel('Frequency'); 
ylabel('Value'); 
 
% 2.6.Demodulation of the NBFM signal.
diff_FM = diff(FM_t); %Differentiate the time domain modulated signal
 
FM_demod_t = resample(abs(hilbert(diff_FM)), Fs, resampleFrequency); %Resample to get the demodulated time domain signal (bring it back to Fs instead of 5Fc)
 
t = linspace(0, length(FM_demod_t)/Fs, length(FM_demod_t)); %Time range to plot the demodulated signal in time domain
figure(fig_i);
fig_i = fig_i + 1;
plot(t, FM_demod_t)
title('Demodulated signal at time domain');
xlabel('Time'); 
ylabel('Value'); 
sound(FM_demod_t, Fs); %Sound the demodulated signal
