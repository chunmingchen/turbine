function plot_fft(time, data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT
Y = fft(data(1,:));

% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
L = length(time);             % Length of signal
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = L;            % Sampling frequency
% f = Fs*(0:(L/2))/L/(L/3600);

T=3600
f = (0:1:(L/2))/(L/T);
% t = 1./f(2:end);

figure
plot(f, P1)

title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (freq/Revolution)')
ylabel('|P1(f)|')


%removing offset in the pressure data
data = data - mean(data);

%figure
% spectrogram(data(1,:),256, 255, 3600, 'yaxis')

end