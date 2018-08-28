function [Px,w] = per_smooth(x,wind,M)
% EE-658 Advanced Digital Signal Processing
% Blackman-Tukey Method
% x: Input signal
% wind: The window sequence
% M: Window's half length

x=x(:);
r=1/M*xcorr(x(1:M),x(1:M));
M=2*M-1;
r=r.*wind;
[Pxx,w]=freqz(r,1);
Px=abs(Pxx);





