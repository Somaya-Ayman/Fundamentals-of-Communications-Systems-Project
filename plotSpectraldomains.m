function plotSpectraldomains = plotSpectraldomains (lineCodedStream, frequency, ylimit)
% fourier transform and normalize the signal
NN = length(lineCodedStream);
fftbitSream = fftshift(fft(lineCodedStream));
% plot the spectral domain
figure
plot(frequency,abs(fftbitSream.^2)/NN);
grid on;
xlabel("Frequency");
ylabel("Power Spectral density" );
ylim([0 ylimit]);
xlim([-0.1 0.1]);
end
