% returns the cutoff time tc at which the trailing pulse 
% envelope falls below tpe dB with respect to the peak envelope amplitude.
tc = gauspuls('cutoff',50e3,0.6,[],-40); 

% define time array needed
t = -tc : 1e-7 : tc; 

% generate pulse, hat 
[yi,yq,ye] = gauspuls(t,50e3,0.6); 

plot(t,yi,t,yq,t,ye)
legend('Inphase','Quadrature','Envelope')
%% Pulse Generation

% Define variables
f0 = 1e6;       %[Hz]
bw = 0.25*f0;   %[Hz]
% Minimum and maximum supportable frequencies
fmin = 0; 
fmax = 2*f0;
w = 2* pi *linspace(fmin,fmax,1000); 

% Define phase
phi = pi/2;

% Plot Gaussian spectrum
F = gaus(w, 2*pi*f0, bw)* exp(1i*phi);
plot(w, F);

% Fourier transform and plot transformed pulse
f = ifft(F);
dt = 1/fmax; 
t = (0:length(f)-1)*dt;

G = gaus(w, 2*pi*f0, bw);
plot(ifftshift(t), burst(G, 0), ifftshift(t), burst(G,pi/2), ifftshift(t), burst(G,pi/4))
legend('0','pi/2','pi/4')
title('Phase shift and Pulse shape')
xlabel('Time t /s')
%% Creates a series of burst pulses with different phases
burst(G, [0,1,2])

%% Function Definitions
function y = gaus(x, mu, sd)
% GAUS generates a normalised Gaussian signal
% x is the sampling domain
% mu is the location of the centre
% sd is the standard deviation 
    y = 1/(2*pi*sd)* exp(-(x-mu).^2/(2*sd^2));
end


function f = burst(G, phi)
% BURST creates a pulse with spectrum G in the Fourier domain
% but is entirely phased shifted by phase phi
    F = G.* exp(1i*phi).';
    f = ifft(F); 
end


