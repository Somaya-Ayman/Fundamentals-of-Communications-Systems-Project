%Transmitter

% generate stream of 100 random bits
bits = randi([0 1],1,100);
% make the zeros into -1 for polar non return to zero line coding
for index = 1:length(bits)
  if bits(index) == 0
    bits(index)= -1;
  endif
end
% extend the duration of each bit (each bit will be represented by 100 values)
delay = ones(100,1);
pulses = bits.* delay;
pulses = reshape(pulses,1,[]);

%time domain
N = length(pulses);
ts = 1e-11;
T = ts * N;
tb = 100* ts;
t = 0:ts:T-ts;

%frequency domain
df = 1/T;
fs = 1/ts;
Rb = 1/tb;
%N is even
f = -(0.5 * fs):df: (0.5 * fs -df);

% fourier transform and normalize the signal
PULSES = fftshift(fft(pulses))/N;
% plot the spectral domain
PULSES = PULSES .* PULSES;
figure
plot(f,abs(PULSES)); grid on; xlabel("Frequency"); ylabel("Normalized Power Spectral density"); xlim([-1e10 1e10]);

%BPSK modulation
fc = 1e9;
carrier = cos(2 * pi * fc * t);
modPulses =  pulses .* carrier;
%plot modulated BPSK with time
figure
plot(t,modPulses); grid on; xlabel("Time"); ylabel("Modulated BPSK signal"); xlim([0 10e-9]); ylim([-1.1 1.1]);

%Plot spectrum of the modulated BPSK signal
MODPULSES = fftshift(fft(modPulses))/N ;
MODPULSES = MODPULSES .* MODPULSES;
figure
plot(f,abs(MODPULSES)); grid on; xlabel("Frequency"); ylabel("Spectrum of modulated BPSK signal");xlim([-1e10 1e10]);


%Receiver
carrier2 = cos(2 * pi * fc * t);
modPulses2 =  modPulses .* carrier2;
MODPULSES2 = fftshift(fft(modPulses2))/N;
MODPULSES2 = MODPULSES2 .* MODPULSES2;

%low pass filter
h = abs(f)<1e9+50;
MODPULSES2 = h.* MODPULSES2;
MODPULSES2 = sqrt(MODPULSES2);
pulses2 = real(ifft(fftshift(MODPULSES2))* N);

average = ones(1,100);
for index = 1:length(average)
  average(index) = pulses(100 * index -50);     % get mid point of each pulse
end
%Decision device
decision = sign(average);

%Compare output of decision level with the generated stream of bits
errorBits = 0;
for index = 1:length(bits)
   if bits(index) != decision(index)
      errorBits++;
  endif
end
bitErrorRate = errorBits/length(bits);
