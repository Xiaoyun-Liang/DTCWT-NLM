function res=epsi(SNR)

% Pierrick Coupe - pierrick.coupe@gmail.com                                                                        
% Brain Imaging Center, Montreal Neurological Institute.                     
% Mc Gill University                                                         
%                                                                            
% Copyright (C) 2008 Pierrick Coupe 

% Based on Koay estimation of truth SNR from Magnitude SNR.

if (SNR > 37) res =1;

else
    
res = 2 + SNR^2 - pi/8 * exp(-(SNR^2/2))*((2+SNR^2)*besseli(0,(SNR^2/4)) + SNR^2*besseli(1,(SNR^2/4)))^2;

end