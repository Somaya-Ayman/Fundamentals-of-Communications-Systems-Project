%part one
%transmitter

bitsStream = randi([0 1],1,10e3);

%unipolar NRZ line coding
positiveVoltage = 1.2;
UniPolarNRZbitSream = positiveVoltage*bitsStream;

%extend the duration of each bit (each bit will be represented by 100 values)
delay = ones(100,1);
UniPolarNRZbitSream = UniPolarNRZbitSream.* delay;
UniPolarNRZbitSream = reshape(UniPolarNRZbitSream,1,[]);

%time domain
N = length(UniPolarNRZbitSream);
tb = 100;
ts = tb/100;
T = ts * N;
t = 0:ts:((N-1)*ts);
n = length(bitsStream);

%frequency domain
df = 1/T;
Rb = 1/tb;
fs = 1/ts;
if(rem(N,2)==0) %% Even
    f = - (0.5*fs) : df : (0.5*fs-df) ;
else %% Odd
    f = - (0.5*fs-0.5*df) : df : (0.5*fs-0.5*df)
end


%plot the UniPolarNRZbitSream
figure;
plot(t,UniPolarNRZbitSream);
grid on; xlabel("Time"); ylabel("UniPolar NRZ Stream"); xlim([0 1000]);ylim([-2 2])
%eyediagram
pkg load communications
eyediagram(UniPolarNRZbitSream, 300); ylabel("UniPolarNRZbitSream eyediagram")

% Compute the Fourier transform

plotSpectraldomains(UniPolarNRZbitSream,f, 350); ylabel("unipolar nrz power spectrum")


%part one
%receiver

% Line encoding using Polar NRZ

polar_nrz = positiveVoltage*(bitsStream*2-1); % Convert 0's to -1.2's and 1's to +1.2's


%extend the duration of each bit (each bit will be represented by 100 values)
delay = ones(100,1);
polar_nrz = polar_nrz.* delay;
polar_nrz = reshape(polar_nrz,1,[]);

%plot the polar_nrz
figure;
plot(t,polar_nrz);
grid on; xlabel("Time"); ylabel("polar NRZ Stream"); xlim([0 1000]);ylim([-2 2])


%eye diagram
eyediagram(polar_nrz, 500); ylabel("polar nrz eyediagram")


% Compute the Fourier transform

% fourier transform and normalize the signal
plotSpectraldomains(polar_nrz,f,1300); ylabel("polar nrz power spectrum")

% Line encoding using Uni-polar RZ
unipolar_rz = positiveVoltage*bitsStream;

%extend the duration of each bit (each bit will be represented by 100 values)
delay = ones(100,1);
unipolar_rz = unipolar_rz.* delay;
unipolar_rz = reshape(unipolar_rz,1,[]);

% Line encoding using Uni-polar RZ
for i = 1:length(unipolar_rz)
    if unipolar_rz(i) == positiveVoltage
      i = i +50;
      j = 0;
      unipolar_rz(i) = 0;
      j++;
      if j == 50
        continue;
      end

    else
      unipolar_rz(i) = 0;
end
end


%plot the unipolar_rz
figure;
plot(t,unipolar_rz);
grid on; xlabel("Time"); ylabel("unipolar RZ Stream"); xlim([0 1000]);ylim([-2 2])

%eye diagram unipolar rz
eyediagram(unipolar_rz, 500); ylabel("unipolar rz eyediagram")


% Compute the Fourier transform

% fourier transform and normalize the signal
plotSpectraldomains(unipolar_rz,f,100); ylabel("unipolar rz power spectrum")



% Line encoding using Bipolar RZ

prev_polarity = 1;
for i = 1:length(bitsStream)
    if bitsStream(i) == 1
        BipolarbitsStream(i) = bitsStream(i)*prev_polarity;
        prev_polarity = -prev_polarity;
    else
        BipolarbitsStream(i) = 0;
    end
end

%extend the duration of each bit (each bit will be represented by 100 values)
delay = ones(100,1);
bipolar_rz = BipolarbitsStream*positiveVoltage;
bipolar_rz = bipolar_rz.* delay;
bipolar_rz = reshape(bipolar_rz,1,[]);

% Line encoding using Bipolar RZ
for i = 1:length(bipolar_rz)
    if bipolar_rz(i) != 0
      i = i +50;
       j = 0;
      bipolar_rz(i) = 0;
      j++;
      if j == 50
        continue;
      end
    else
      bipolar_rz(i) = 0;
end
end


%plot the bipolar_rz
figure;
plot(t,bipolar_rz);
grid on; xlabel("Time"); ylabel("bipolar RZ Stream"); xlim([0 1000]);ylim([-2 2])

%eye diagram bipolar_rz
eyediagram(bipolar_rz, 500); ylabel("bipolar RZ eyediagram")

% Compute the Fourier transform bipolar rz

% fourier transform and normalize the signal
plotSpectraldomains(bipolar_rz,f,300); ylabel("bipolar rz power spectrum")


% Line encoding using Manchester
Manchester = positiveVoltage*(bitsStream*2-1); % Convert 0's to -1.2's and 1's to +1.2's

%extend the duration of each bit (each bit will be represented by 100 values)
delay = ones(100,1);
Manchester = Manchester.* delay;
Manchester = reshape(Manchester,1,[]);

% Line encoding using Manchester
for i = 1:length(bitsStream)
    if bitsStream(i) == 0
        Manchester((i-1)*100+1:50+(i-1)*100) = -1.2;
        Manchester(51+(i-1)*100:i*100) = 1.2;
     else
        Manchester((i-1)*100+1:50+(i-1)*100) = 1.2;
        Manchester(51+(i-1)*100:i*100) = -1.2;
    end
end

%plot the Manchester
figure;
plot(t,Manchester(1:N));
grid on; xlabel("Time"); ylabel("Manchester Stream"); xlim([0 1000]);ylim([-2 2])

%eye diagram manchester
eyediagram(Manchester, 500); ylabel("Manchester eyediagram")

% Compute the Fourier transform of manchester

% fourier transform and normalize the signal
plotSpectraldomains(Manchester,f,600); ylabel("Manchester power spectrum")



% Receiver code for decoding the received signal and calculating BER
%before adding the noise
rx_unipolar_nrz = UniPolarNRZbitSream;
rx_polar_nrz = polar_nrz;
rx_unipolar_rz = unipolar_rz;
rx_bipolar_rz = bipolar_rz;
rx_manchester = Manchester;


%tx bitstream extended
delay = ones(100,1);
extended_bit_stream = bitsStream.* delay;
extended_bit_stream = reshape(extended_bit_stream,1,[]);


%decision making

% Decoding the received signal using Uni-polar RZ
rx_bitsAfterDecision_unipolar_rz = zeros(1,length(rx_unipolar_rz));
i = 1;
for i = 1:100:length(rx_unipolar_rz)
    if rx_unipolar_rz(i) > 0.6
       rx_bitsAfterDecision_unipolar_rz(i:i+99) = 1;
    else
       rx_bitsAfterDecision_unipolar_rz(i:i+99) = 0;
    end
    i = i + 1;
end

% Decoding the received signal using Uni-polar NRZ
rx_bitsAfterDecision_unipolar_nrz = zeros(1,length(rx_unipolar_nrz));
i = 1;
for i = 1:length(rx_unipolar_nrz)
    if rx_unipolar_nrz(i) > 0.6
       rx_bitsAfterDecision_unipolar_nrz(i) = 1;
    else
       rx_bitsAfterDecision_unipolar_nrz(i) = 0;
    end
    i = i + 1;
end

% Decoding the received signal using polar NRZ
rx_bitsAfterDecision_polar_nrz = zeros(1,length(rx_polar_nrz));
i = 1;
for i = 1:100:length(rx_polar_nrz)
    if rx_polar_nrz(i) > 0
       rx_bitsAfterDecision_polar_nrz(i:i+99) = 1;
    else
       rx_bitsAfterDecision_polar_nrz(i:i+99) = 0;
    end
    i = i + 1;
end

% Decoding the received signal using bi-polar RZ
rx_bitsAfterDecision_bipolar_rz = zeros(1,length(rx_bipolar_rz));

i = 1;
for i = 1:100:length(rx_bipolar_rz)
    if rx_bipolar_rz(i) > 0.6 || rx_bipolar_rz(i) < -0.6
       rx_bitsAfterDecision_bipolar_rz(i:i+99) = 1;
    else
       rx_bitsAfterDecision_bipolar_rz(i:i+99) = 0;
    end
    i = i + 1;
end

% Decoding the received signal using manchester
rx_bitsAfterDecision_manchester = zeros(1,length(rx_manchester));
i = 1;
for i = 1:100:length(rx_manchester)
    if rx_manchester(i) >  0
       rx_bitsAfterDecision_manchester(i:i+99) = 1;
    else
       rx_bitsAfterDecision_manchester(i:i+99) = 0;
    end
    i = i + 1;
end

% Calculating BER for each line coding technique
ber_unipolar_nrz = sum(rx_bitsAfterDecision_unipolar_nrz~=extended_bit_stream)/length(extended_bit_stream);
ber_polar_nrz = sum(rx_bitsAfterDecision_polar_nrz~=extended_bit_stream)/length(extended_bit_stream);
ber_unipolar_rz = sum(rx_bitsAfterDecision_unipolar_rz~=extended_bit_stream)/length(extended_bit_stream);
ber_bipolar_rz = sum(rx_bitsAfterDecision_bipolar_rz~=extended_bit_stream)/length(extended_bit_stream);
ber_manchester = sum(rx_bitsAfterDecision_manchester~=extended_bit_stream)/length(extended_bit_stream);
% Displaying the results
disp(['BER for Uni-polar NRZ = ' num2str(ber_unipolar_rz)]);
disp(['BER for Polar NRZ = ' num2str(ber_polar_nrz)]);
disp(['BER for Uni-polar RZ = ' num2str(ber_unipolar_rz)]);
disp(['BER for Bipolar RZ = ' num2str(ber_bipolar_rz)]);
disp(['BER for Manchester = ' num2str(ber_manchester)]);

%define sigma
sigma = linspace(0, 1.2, 10); % Generate 10 equally spaced values between 0 and 1.2


%decision making

% Decoding the received signal using Uni-polar RZ
rx_bitsAfterDecision_unipolar_rz = zeros(1,length(rx_unipolar_rz));
ber_unipolar_rz = zeros(1, 10);
for j = 1:length(sigma)
rx_unipolar_rz =  sigma(j).* randn(1,length(t)) + unipolar_rz;
i = 1;
for i = 1:100:length(rx_unipolar_rz)
    if rx_unipolar_rz(i) > 0.6
       rx_bitsAfterDecision_unipolar_rz(i:i+99) = 1;
    else
       rx_bitsAfterDecision_unipolar_rz(i:i+99) = 0;
    end
    i = i + 1;
end
ber_unipolar_rz(j) = sum(rx_bitsAfterDecision_unipolar_rz~=extended_bit_stream)/length(extended_bit_stream);
end
%plot semilogy
figure;
semilogy(sigma, ber_unipolar_rz);xlabel('Noise RMS value');
ylabel('Bit error rate (BER) unipolar rz');

% Decoding the received signal using Uni-polar NRZ
rx_bitsAfterDecision_unipolar_nrz = zeros(1,length(rx_unipolar_nrz));
ber_unipolar_nrz = zeros(1, 10);
for j = 1:length(sigma)
rx_unipolar_nrz =  sigma(j).* randn(1,length(t)) + UniPolarNRZbitSream;
i = 1;
for i = 1:length(rx_unipolar_nrz)
    if rx_unipolar_nrz(i) > 0.6
       rx_bitsAfterDecision_unipolar_nrz(i) = 1;
    else
       rx_bitsAfterDecision_unipolar_nrz(i) = 0;
    end
    i = i + 1;
end
ber_unipolar_nrz(j) = sum(rx_bitsAfterDecision_unipolar_nrz~=extended_bit_stream)/length(extended_bit_stream);
end
%plot semilogy
figure;
semilogy(sigma, ber_unipolar_nrz);xlabel('Noise RMS value');
ylabel('Bit error rate (BER) unipolar nrz');

% Decoding the received signal using polar NRZ
rx_bitsAfterDecision_polar_nrz = zeros(1,length(rx_polar_nrz));
ber_polar_nrz = zeros(1, 10);
for j = 1:length(sigma)
rx_polar_nrz =  sigma(j).* randn(1,length(t)) + polar_nrz;
i = 1;
for i = 1:100:length(rx_polar_nrz)
    if rx_polar_nrz(i) > 0
       rx_bitsAfterDecision_polar_nrz(i:i+99) = 1;
    else
       rx_bitsAfterDecision_polar_nrz(i:i+99) = 0;
    end
    i = i + 1;
    end
    ber_polar_nrz(j) = sum(rx_bitsAfterDecision_polar_nrz~=extended_bit_stream)/length(extended_bit_stream);
end
%plot semilogy
figure;
semilogy(sigma, ber_polar_nrz);xlabel('Noise RMS value');
ylabel('Bit error rate (BER) polar nrz');


% Decoding the received signal using bi-polar RZ
rx_bitsAfterDecision_bipolar_rz = zeros(1,length(rx_bipolar_rz));
%save values for rx bits with noise and rx bits after decision
rxBiPolar = zeros(10, length(t));
rxBiPolarDecision = zeros(10, length(t));
ber_bipolar_rz = zeros(1, 10);
for j = 1:length(sigma)
rx_bipolar_rz =  sigma(j).* randn(1,length(t)) + bipolar_rz;
rxBiPolar(j,:) = rx_bipolar_rz;
i = 1;
for i = 1:100:length(rx_bipolar_rz)
    if rx_bipolar_rz(i) > 0.6 || rx_bipolar_rz(i) < -0.6
       rx_bitsAfterDecision_bipolar_rz(i:i+99) = 1;
    else
       rx_bitsAfterDecision_bipolar_rz(i:i+99) = 0;
    end
    i = i + 1;
  end
  rxBiPolarDecision(j,:) = rx_bitsAfterDecision_bipolar_rz;
  ber_bipolar_rz(j) = sum(rx_bitsAfterDecision_bipolar_rz~=extended_bit_stream)/length(extended_bit_stream);
end
%plot semilogy
figure;
semilogy(sigma, ber_bipolar_rz);xlabel('Noise RMS value');
ylabel('Bit error rate (BER) bipolar rz');

% Decoding the received signal using manchester
rx_bitsAfterDecision_manchester = zeros(1,length(rx_manchester));
ber_manchester = zeros(1, 10);
for j = 1:length(sigma)
rx_manchester =  sigma(j).* randn(1,length(t)) + Manchester;
i = 1;
for i = 1:100:length(rx_manchester)
    if rx_manchester(i) >  0
       rx_bitsAfterDecision_manchester(i:i+99) = 1;
    else
       rx_bitsAfterDecision_manchester(i:i+99) = 0;
    end
    i = i + 1;
  end
  ber_manchester(j) = sum(rx_bitsAfterDecision_manchester~=extended_bit_stream)/length(extended_bit_stream);

end
%plot semilogy
figure;
semilogy(sigma, ber_manchester);xlabel('Noise RMS value');
ylabel('Bit error rate (BER) manchester');

% Calculating BER for each line coding technique
% Displaying the results
disp(['BER for Uni-polar NRZ with noise = ' num2str(ber_unipolar_rz)]);
disp(['BER for Polar NRZ with noise = ' num2str(ber_polar_nrz)]);
disp(['BER for Uni-polar RZ with noise = ' num2str(ber_unipolar_rz)]);
disp(['BER for Bipolar RZ with noise = ' num2str(ber_bipolar_rz)]);
disp(['BER for Manchester with noise = ' num2str(ber_manchester)]);
%bonus
% Error detection

detected_bits = zeros(1, 10);
for j = 1:length(sigma)
    previousPolarity = -1;
    for i = 50:100:length(t) %start from mid of bit
        if rxBiPolar(j,i) > 0.6 || rxBiPolar(j,i) < -0.6 %check if consecutive ones have same polarity
            if sign(rxBiPolar(j,i)) == previousPolarity
                detected_bits(j) = detected_bits(j) + 1;
            end
            previousPolarity = sign(rxBiPolar(j,i));
        end
    end
    disp(['Number of detected errors for sigma ' num2str(sigma(j)) ' is ' num2str(detected_bits(j))]);
end
%plot semilogy
figure;
semilogy(sigma, detected_bits/length(extended_bit_stream));xlabel('Noise RMS value');
ylabel('Bit error rate (BER) bipolar rz bonus part');


