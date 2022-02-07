tc = gauspuls('cutoff',50e3,0.6,[],-40); 
t = -tc : 1e-7 : tc; 
[yi,yq,ye] = gauspuls(t,50e3,0.6); 

plot(t,yi,t,yq,t,ye)
legend('Inphase','Quadrature','Envelope')