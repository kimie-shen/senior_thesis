Fs = 44100;         % Sampling frequency
T = 1/Fs;           % Sampling period
t = -0.5:T:0.5;     % Time vector
L = length(t);      % Signal length

X = 1/(0.4*sqrt(2*pi))*(exp(-t.^2/(2*(0.1*1e-3)^2)));

plot(t,X)
title("Gaussian Pulse in Time Domain")
xlabel("Time (t)")
ylabel("X(t)")
axis([-1e-3 1e-3 0 1.1]) 
n = 2^nextpow2(L);
Y = fft(X,n);
f = Fs*(0:(n/2))/n;
P = abs(Y/sqrt(n)).^2;

plot(f,P(1:n/2+1)) 
title("Gaussian Pulse in Frequency Domain")
xlabel("f (Hz)")
ylabel("|P(f)|")