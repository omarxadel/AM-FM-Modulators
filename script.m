% VARIABLES
fig_i = 1; %figures counter
% 1.Amplitude Modulation
% 1.1.Reading the audio input.
[signal,Fs] = audioread('eric.wav'); 

FrequencyDomainSignal = fftshift(fft(signal)); % input signal in freq domain
f = linspace(-Fs/2,Fs/2,length(FrequencyDomainSignal));
figure(fig_i);
fig_i = fig_i + 1;
plot(f, abs(FrequencyDomainSignal));
title('Spectrum of the input signal');

% 1.2.Filtering frequencies > 4KHz.
BW = 4000; % specified bandwidth
filter = ones(length(FrequencyDomainSignal), 1);
   % creating the ideal filter
    for i = 1:length(FrequencyDomainSignal)
        if f(i)<-BW || f(i)>BW
            filter(i) = 0;
        end
    end
% filtering the signal
FilteredSignal_f = filter.*FrequencyDomainSignal;

figure(fig_i);
fig_i = fig_i + 1;
f_filtered = linspace(-Fs/2,Fs/2,length(real(FilteredSignal_f)));
plot(f_filtered, abs(FilteredSignal_f));
title('Filtered signal in frequency domain');

% 1.3.Filtered signal in time domain.
FilteredSignal_t = real(ifft(ifftshift(FilteredSignal_f)));

figure(fig_i);
fig_i = fig_i + 1;
t = linspace(0,length(FilteredSignal_t)/Fs, length(FilteredSignal_t));
plot(t, real(FilteredSignal_t));
title('Filtered signal in time domain');

% 1.4.Playing the filtered audio signal (sound).
sound(FilteredSignal_t, Fs);
pause(8);

% 1.5.Modulating the signal.
Fc = 100000; % specified carrier frequency
resampleFrequency = 5*Fc;
resampledSignal = resample(FilteredSignal_t, resampleFrequency,Fs);

carrierTime = linspace(0, length(resampledSignal)/resampleFrequency, length(resampledSignal));
carrier_t = cos(2*pi*Fc*carrierTime).'; % creating carrier
carrier_f = fftshift(fft(carrier_t));

DSBSC_t = carrier_t.*resampledSignal; % s(t) = m(t)*c(t)
DSBSC_f = fftshift(fft(DSBSC_t));

figure(fig_i);
fig_i = fig_i + 1;
f_DSBSC=linspace(-resampleFrequency/2,resampleFrequency/2,length(real(DSBSC_f)));
plot(f_DSBSC, abs(DSBSC_f));
title('DSBSC modulated signal in frequency domain');

signalMax = max(abs(resampledSignal));
A = 2*signalMax; %because modulation index =0.5

DSBTC_t = (A+resampledSignal).*carrier_t; % s(t) = (Ac + m(t))*c(t)
DSBTC_f = fftshift(fft(DSBTC_t));

figure(fig_i);
fig_i = fig_i + 1;
f_DSBTC = linspace(-resampleFrequency/2, resampleFrequency/2, length(real(DSBTC_f)));
plot(f_DSBTC, abs(DSBTC_f));
title('DSBTC modulated signal in frequency domain');

% 1.6.Envelope detector demodulation.
[DSBSC_t_demod, DSBSC_f_demod] = envelope_detection(DSBSC_t, Fs, resampleFrequency);
[DSBTC_t_demod, DSBTC_f_demod] = envelope_detection(DSBTC_t, Fs, resampleFrequency);

% 1.7.Sketching and playing the recieved signal.
figure(fig_i);
fig_i = fig_i + 1;
subplot(2,1,1);
t_mod = linspace(0, length(DSBSC_t_demod)/Fs, length(DSBSC_t_demod));
plot(t_mod, real(DSBSC_t_demod));
title('DSBSC demodulated signal in time domain using envelope detector');

subplot(2,1,2);
f_DSBSC_mod = linspace(-Fs/2, Fs/2, length(real(DSBSC_f_demod)));
plot(f_DSBSC_mod, abs(DSBSC_f_demod));
title('DSBSC demodulated signal in frequency domain using envelope detector');

sound(DSBSC_t_demod,Fs);
pause(8);

figure(fig_i);
fig_i = fig_i + 1;
subplot(2,1,1);
t_modTc = linspace(0,length(DSBTC_t_demod)/Fs, length(DSBTC_t_demod));
plot(t_modTc, real(DSBTC_t_demod));
title('DSBTC demodulated signal in time domain using envelope detector');

subplot(2,1,2);
f_DSBTC_mod = linspace(-Fs/2,Fs/2,length(real(DSBTC_f_demod)));
plot(f_DSBTC_mod,abs(DSBTC_f_demod));
title('DSBTC demodulated signal in frequency domain using envelope detector');

sound(DSBTC_t_demod,Fs); 
pause(8);

% 1.8.Adding noise then demodulating the DSB-TC signal.
snr = 0;
DSBTC_t_noise = awgn(DSBTC_t, snr);
[DSBTC_t_demod, ~] = envelope_detection(DSBTC_t_noise, Fs, resampleFrequency);
figure(fig_i);
fig_i = fig_i + 1;
t_modTc = linspace(0,length(DSBTC_t_demod)/Fs, length(DSBTC_t_demod));
plot(t_modTc, real(DSBTC_t_demod));
title({'DSBTC demodulated signal in time domain', 'using envelope detector with snr=0'});

sound(DSBTC_t_demod, Fs);
pause(8);

snr = 10;
DSBTC_t_noise = awgn(DSBTC_t, snr);
[DSBTC_t_demod, ~] = envelope_detection(DSBTC_t_noise, Fs, resampleFrequency);
figure(fig_i);
fig_i = fig_i + 1;
t_modTc = linspace(0,length(DSBTC_t_demod)/Fs, length(DSBTC_t_demod));
plot(t_modTc, real(DSBTC_t_demod));
title({'DSBTC demodulated signal in time domain', 'using envelope detector with snr=10'});

sound(DSBTC_t_demod, Fs);
pause(8);

snr = 30;
DSBTC_t_noise = awgn(DSBTC_t, snr);
[DSBTC_t_demod, ~] = envelope_detection(DSBTC_t_noise, Fs, resampleFrequency);
figure(fig_i);
fig_i = fig_i + 1;
t_modTc = linspace(0,length(DSBTC_t_demod)/Fs, length(DSBTC_t_demod));
plot(t_modTc, real(DSBTC_t_demod));
title({'DSBTC demodulated signal in time domain', 'using envelope detector with snr=30'});

sound(DSBTC_t_demod, Fs);
pause(8);

% 1.9.Coherent detection of the DSB-SC;
[coh_op_t, coh_op_f] = coherent_detection(DSBSC_t, carrier_t, Fs, resampleFrequency, BW);

t_modSc = linspace(0,length(coh_op_t)/Fs, length(coh_op_t));
figure(fig_i);
fig_i = fig_i + 1;
subplot(2,1,1);
plot(t_modSc, real(coh_op_t));
title('DSBSC demodulated signal in time domain using coherent detection');
sound(real(coh_op_t), Fs);
pause(8);

f_DSBSC_mod = linspace(-Fs/2,Fs/2,length(real(coh_op_f)));
subplot(2,1,2);
plot(f_DSBSC_mod, abs(coh_op_f));
title('DSBSC demodulated signal in frequency domain using coherent detection');

snr = 0;
DSBSC_t_noise = awgn(DSBSC_t, snr);
[coh_op_t, coh_op_f] = coherent_detection(DSBSC_t_noise, carrier_t, Fs, resampleFrequency, BW);

t_modSc = linspace(0,length(coh_op_t)/Fs, length(coh_op_t));
figure(fig_i);
fig_i = fig_i + 1;
subplot(2,1,1);
plot(t_modSc, real(coh_op_t));
title({'DSBSC demodulated signal in time domain' ,'using coherent detection with snr = 0'});

f_DSBSC_mod = linspace(-Fs/2,Fs/2,length(real(coh_op_f)));
subplot(2,1,2);
plot(f_DSBSC_mod, abs(coh_op_f));
title({'DSBSC demodulated signal in frequency domain',' using coherent detection with snr = 0'});
sound(real(coh_op_t), Fs);
pause(8);

snr = 10;
DSBSC_t_noise = awgn(DSBSC_t, snr);
[coh_op_t, coh_op_f] = coherent_detection(DSBSC_t_noise, carrier_t, Fs, resampleFrequency, BW);

t_modSc = linspace(0,length(coh_op_t)/Fs, length(coh_op_t));
figure(fig_i);
fig_i = fig_i + 1;
subplot(2,1,1);
plot(t_modSc, real(coh_op_t));
title({'DSBSC demodulated signal in time domain',' using coherent detection with snr = 10'});

f_DSBSC_mod = linspace(-Fs/2,Fs/2,length(real(coh_op_f)));
subplot(2,1,2);
plot(f_DSBSC_mod, abs(coh_op_f));
title({'DSBSC demodulated signal in frequency domain',' using coherent detection with snr = 10'});
sound(real(coh_op_t), Fs);
pause(8);

snr = 30;
DSBSC_t_noise = awgn(DSBSC_t, snr);
[coh_op_t, coh_op_f] = coherent_detection(DSBSC_t_noise, carrier_t, Fs, resampleFrequency, BW);

t_modSc = linspace(0,length(coh_op_t)/Fs, length(coh_op_t));
figure(fig_i);
fig_i = fig_i + 1;
subplot(2,1,1);
plot(t_modSc, real(coh_op_t));
title({'DSBSC demodulated signal in time domain',' using coherent detection with snr = 30'});

f_DSBSC_mod = linspace(-Fs/2,Fs/2,length(real(coh_op_f)));
subplot(2,1,2);
plot(f_DSBSC_mod, abs(coh_op_f));
title({'DSBSC demodulated signal in frequency domain',' using coherent detection with snr = 30'});
sound(real(coh_op_t), Fs);
pause(8);

% 1.10.Coherent detection of the DSB-SC with frequency error.
Fc_error = 100100; % Frequency error of .1 KHz
carrier_error = cos(2*pi*Fc_error*carrierTime).';

[coh_op_t, coh_op_f] = coherent_detection(DSBSC_t, carrier_error, Fs, resampleFrequency, BW);

t_modSc = linspace(0,length(coh_op_t)/Fs, length(coh_op_t));
figure(fig_i);
fig_i = fig_i + 1;
subplot(2,1,1);
plot(t_modSc, real(coh_op_t));
title({'DSBSC demodulated signal in time domain', 'using coherent detection with freq error'});
sound(real(coh_op_t), Fs);
pause(8);

f_DSBSC_mod = linspace(-Fs/2,Fs/2,length(real(coh_op_f)));
subplot(2,1,2);
plot(f_DSBSC_mod, abs(coh_op_f));
title({'DSBSC demodulated signal in frequency domain',' using coherent detection with freq error'});

% 1.11.Coherent detection of the DSB-SC with phase error.
carrier_error = cos((2*pi*Fc*carrierTime) + 20).'; % Phase error 20°

[coh_op_t, coh_op_f] = coherent_detection(DSBSC_t, carrier_error, Fs, resampleFrequency, BW);

t_modSc = linspace(0,length(coh_op_t)/Fs, length(coh_op_t));
figure(fig_i);
fig_i = fig_i + 1;
subplot(2,1,1);
plot(t_modSc, real(coh_op_t));
title({'DSBSC demodulated signal in time domain using coherent detection',' with phase error'});
sound(real(coh_op_t), Fs);
pause(8);

f_DSBSC_mod = linspace(-Fs/2,Fs/2,length(real(coh_op_f)));
subplot(2,1,2);
plot(f_DSBSC_mod, abs(coh_op_f));
title({'DSBSC demodulated signal in frequency domain using coherent detection',' with phase error'});

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
