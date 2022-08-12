clear;clc;
close all

sgn = @(in) sign(in+eps);

K = 3;      % Sparsity
N = 5000;   % Signal dimension  
fs = 5000;  % Sampling frequency  
ts = 1/fs;  % Sampling interval  
Ts = 1:N;   % Sampling sequence 
original_s = 0.9*cos(100*pi*Ts*ts)+0.2*cos(150*pi*Ts*ts)+cos(200*pi*Ts*ts);  %real signal
x_bit = sgn(original_s);                    % 1-bit real signal
Psi = fft(eye(N,N))/sqrt(N);                % Fourier positive transform matrix  
T = Psi';                                   % Recovery matrix (measurement matrix * orthogonal inverse transformation matrix) 
spect= sub_BIHT_l1(T,x_bit',2*K);           % Get the normalized frequency spectrum

figure 
plot(abs(spect))
xlabel('Frequency(Hz)')
ylabel('Amplitude')
title('The normalized frequency spectrum')

spect(ceil(length(spect)/2):end) = 0;        % Eliminating the negative frequencies
x_imag_bit1 = sgn(imag(Psi'*spect));         % 1-bit Hilbert Transform

figure
plot(Ts*ts,x_imag_bit1,'b:');
xlabel('Time(s)')
ylabel('Amplitude')
legend('$\mathcal{H}_b\textnormal{(sgn}(u(t)))$')

x_imag = imag(hilbert(original_s));                
x_imag_bit2 = sgn(x_imag);                          %  metric           
x_imag_bit3 = sgn(imag(hilbert(sgn(original_s))));  %  handled as square wave
x_imag_square = imag(hilbert(sgn(original_s)));     %  handled as square wave

figure
plot(Ts*ts,x_imag_bit2,'r-')
xlabel('Time(s)')
ylabel('Amplitude')
legend('$\textnormal{sgn}(\mathcal{H}(u(t)))$')

hd1 = nnz(x_imag_bit1 - x_imag_bit2')/N             % Normalized Hamming distance
hd2 = nnz(x_imag_bit3 - x_imag_bit2)/N
R1 = corrcoef(x_imag_bit2,x_imag_bit1)              % Correlation coefficient
R2 = corrcoef(x_imag_bit2,x_imag_bit3)

figure
plot(Ts*ts,x_imag_bit2,'r-','linewidth',2)
hold on
plot(Ts*ts,x_imag_square,'k--','linewidth',2)
hold on
plot(Ts*ts,x_imag_bit1,':b','linewidth',2)
axis([0 0.05 -3 3]) ;
xlabel('Time(s)')
ylabel('Amplitude')
grid on
legend('$\textnormal{sgn}(\mathcal{H}(u(t)))$','$\mathcal{H}\textnormal{(sgn}(u(t)))$','$\mathcal{H}_b\textnormal{(sgn}(u(t)))$')