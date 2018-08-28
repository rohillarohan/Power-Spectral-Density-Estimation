function Px = MV(x,p,w)
% EE-658 Advanced Digital Signal Processing
% Minimum Variance Spectrum estimation
% x: Input signal
% p: The order of the filters

N=2*length(w);
x=x(:);
r=1/p*xcorr(x(1:p),x(1:p));
M=length(r);
R=toeplitz(r((M+1)/2:end));
[V,D]=eig(R);
d=1./(abs(diag(D))+eps);
W=abs(fft(V,N)).^2;
Px1=fftshift((1+p)./(W*d));
Px=Px1(N/2+1:end)';
% w=-pi:(2*pi)/N:pi-pi/N;