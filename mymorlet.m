function wavelet= mymorlet(M, w, s)
% M : length of the wavelet
% w: Omega0 = 5.0
% s: scaling factor = 1.0
% frequency in Hz : f = 2*s*w*r / M.  r is the sampling rate

%t= -1:(1/length):1;    
x = linspace(-s * 2 * pi, s * 2 * pi, M);
%s = 6/(2*pi*f)
wavelet= exp(1i * w * x) .* exp(-0.5 * x.^2) * pi^(-0.25);
%exp(2*pi*1i*f.*t)  .* exp(-t.^2./(2*s^2));;

end