[signal,Fs] = audioread('eric.wav'); 

FrequencyDomainSignal = fftshift(fft(signal));
f=linspace(-Fs/2,Fs/2,length(FrequencyDomainSignal));
figure(1);
plot(f,abs(FrequencyDomainSignal))
title('Spectrum of the signal')

BW=4000;
 filter = ones(length(FrequencyDomainSignal),1);
   
    for i = 1: length(FrequencyDomainSignal)
        if f(i)<-BW || f(i)>BW
            filter(i)=0;
        end
    end
    
FilteredSignal_f = filter.*FrequencyDomainSignal;
figure();
f_filtered=linspace(-Fs/2,Fs/2,length(real(FilteredSignal_f)));
plot(f_filtered,abs(FilteredSignal_f));
title('filtered signal in frequency doain');

FilteredSignal_t= real(ifft(ifftshift(FilteredSignal_f)));
sound(FilteredSignal_t,Fs);

figure();
t = linspace(0,length(FilteredSignal_t)/Fs, length(FilteredSignal_t));
plot(t,real(FilteredSignal_t));
title('filtered signal in time doain');

Fc=100000;
resampleFrequency=5*Fc;
resampledSignal=resample(FilteredSignal_t,resampleFrequency,Fs);

carrierTime = linspace(0,length(resampledSignal)/resampleFrequency, length(resampledSignal));
carrier_t = cos(2*pi*Fc*carrierTime).';
carrier_f = fftshift(fft(carrier_t));

DSBSC_t = carrier_t.*resampledSignal;
DSBSC_f = fftshift(fft(DSBSC_t));

figure();
f_DSBSC=linspace(-resampleFrequency/2,resampleFrequency/2,length(real(DSBSC_f)));
plot(f_DSBSC,abs(DSBSC_f));
title('DSBSC modulated signal in frequency domain');
    
signalMax= max(abs(resampledSignal));
A = 2*signalMax;
modulationIndex=0.5;

DSBTC_t = (A+resampledSignal).*carrier_t;
DSBTC_f = fftshift(fft(DSBTC_t));

figure();
f_DSBTC=linspace(-resampleFrequency/2,resampleFrequency/2,length(real(DSBTC_f)));
plot(f_DSBTC,abs(DSBTC_f));
title('DSBTC modulated signal in frequency domain');


snr_SC=0;
 if snr_SC == 1
        DSBSC_t = awgn(DSBSC_t, snr_SC);
 end
    envelope_sc=abs(hilbert(DSBSC_t));
    DSBSC_t_demod = resample(envelope_sc,resampleFrequency,Fs); %envelope detector and resample 
    DSBSC_f_demod = fftshift(fft(DSBSC_t_demod));
    
figure();
t_mod = linspace(0,length(DSBSC_t_demod)/Fs, length(DSBSC_t_demod));
plot(t_mod,real(DSBSC_t_demod));
title('DSBSC demodulated signal in time domain using envelope detector');

figure();
f_DSBSC_mod=linspace(-Fs/2,Fs/2,length(real(DSBSC_f_demod)));
plot(f_DSBSC_mod,abs(DSBSC_f_demod));
title('DSBSC demodulated signal in frequency domain using envelope detector');
%sound(DSBSC_t_demod,Fs); %should not be good

    
snr_TC=0;
 if snr_TC == 1
        DSBTC_t = awgn(DSBTC_t, snr_TC);
 end
    envelope_tc=abs(hilbert(DSBTC_t));
    DSBTC_t_demod = resample(envelope_tc,resampleFrequency,Fs); %envelope detector and resample 
    DSBTC_f_demod = fftshift(fft(DSBTC_t_demod));

figure();
t_modTc = linspace(0,length(DSBTC_t_demod)/Fs, length(DSBTC_t_demod));
plot(t_modTc,real(DSBTC_t_demod));
title('DSBTC demodulated signal in time domain using envelope detector');

figure();
f_DSBTC_mod=linspace(-Fs/2,Fs/2,length(real(DSBTC_f_demod)));
plot(f_DSBTC_mod,abs(DSBTC_f_demod));
title('DSBTC demodulated signal in frequency domain using envelope detector');
%sound(DSBTC_t_demod,Fs); %should not be good
