# Fundamentals-of-Communications-Systems-Project
Overview
This project involves the implementation of a line coding system and binary phase shift-keying (BPSK) transmitter and receiver. The project's goal is to design, examine various line coding schemes, and assess how well they work in the face of noise. There are two sections to the project:

Part I: Implementation of a line coding system and evaluation of its performance under different line coding techniques and noise levels.
Part II: Implementation of a BPSK transmitter and receiver and evaluation of its performance.

# Part I: Line Coding System
# Transmitter
Generate stream of random bits (10,000 bit) (This bit stream should be selected to be random, which means that the type of each bit is randomly selected by the program code to be either ‘1’ or ‘0’).
Line code the stream of bits (pulse shape) according to Uni-polar non return to zero (Supply voltages are: +1.2 V and -1.2V).
Plot the corresponding Eye diagram.
Plot the spectral domains of the pulses (square of the Fourier transform).
# Receiver
Design a receiver which consists of a decision device. (The decision device has two inputs: received waveform).
Compare the output of the decision level with the generated stream of bits in the transmitter. The comparison is performed by comparing the value of each received bit with thecorresponding transmitted bit (step 1) and count number of errors. Then calculate bit error rate (BER) = number of error bits/ Total number of bits.
Repeat the previous steps for different line coding (Polar non return to zero, Uni-polar return to zero, Bipolar return to zero and Manchester coding).
Add noise to the received signal (Hint: use n = sigma * randn(1,length(t) ),where t is time vector and sigma is the noise RMS value).
Sweep on the value of sigma (10 values ranges from 0 to the maximum supply voltage) and calculate the corresponding BER for each value of sigma.
Repeat the previous steps for different line coding and plot BER versus sigma for the different line coding in the same figure, where y-axis is in the log scale (Hint: use semilogy).
(Bonus) For the case of Bipolar return to zero, design an error detection circuit. Count the number of detected errors in case of different number of sigma (Use the output of step 8).
# Part II: BPSK Transmitter and Receiver
# Transmitter
Generate stream of random bits (100 bit) (This bit stream should be selected to be random, which means that the type of each bit is randomly selected by the program code to be either ‘1’ or ‘0’).
Line code the stream of bits (pulse shape) according to Polarnon-return to zero (Maximum voltage +1, Minimum voltage -1).
Plot the spectral domains.
Plot the time domain of the modulated BPSK signal (fc = 1 GHz).
Plot the spectrum of the modulated BPSK signal.
# Receiver
Design a receiver which consists of modulator, integrator (simply LPF) and decision device.
Compare the output of decision level with the generated stream of bits in the transmitter. The comparison is performed by comparing the value of each received bit with the corresponding transmitted bit (step 1) and count number of errors. Then calculate bit error rate (BER) = number of error bits/ Total number of bits.
Required Software and Tools
MATLAB/OCTAVE
Signal Processing Toolbox
