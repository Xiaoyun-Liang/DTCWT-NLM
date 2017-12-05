function [sig]=RicianSTD(ima)

% Pierrick Coupe - pierrick.coupe@gmail.com                                  
% Jose V. Manjon - jmanjon@fis.upv.es                                        
% Brain Imaging Center, Montreal Neurological Institute.                     
% Mc Gill University                                                         
%                                                                            
% Copyright (C) 2008 Pierrick Coupe and Jose V. Manjon  

% References

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Rician noise estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Pierrick Coupe, Jose Manjon, Elias Gedamu, Douglas Arnold,
%     Montserrat Robles, and D. Collins. 
%
%     An object-based method for rician noise estimation in mr images.  
%     
%     Medical Image Computing and Computer-Assisted Intervention MICCAI 2009,
%     volume 5762, chapter 73, pages 601-608. Springer Berlin Heidelberg,
%     Berlin, Heidelberg, 2009.


double(ima);
s=size(ima);

% Estimation of the zeros pading size
p(1) = 2^(ceil(log2(s(1))));
p(2) = 2^(ceil(log2(s(2))));
p(3) = 2^(ceil(log2(s(3))));

% Zeros Pading
pad1 = zeros(p(1),p(2),p(3));
pad1(1:s(1),1:s(2),1:s(3)) = ima(:,:,:);

% Wavelet Transform
[af, sf] = farras;
w1 = dwt3D(pad1,1,af);

% Removing region corresponding to zeros pading
tmp = w1{1}{7};
tmp2 = w1{2};
tmp = tmp(1:round((s(1)-1)/2),1:round((s(2)-1)/2),1:round((s(3)-1)/2));
tmp2 = tmp2(1:round((s(1)-1)/2),1:round((s(2)-1)/2),1:round((s(3)-1)/2));

% Detection of the object in the LLL subband
[mu,mask]=kmeans(tmp2,2);
th=mean(mu);
map = (tmp2(:,:,:)>th);

% Detection of the High gradient area in the LLL subband
[PX,PY,PZ] = gradient(tmp2);
GR = sqrt(PX.^2 + PY.^2 + PZ.^2);
m = median(GR(map));
map2 = (GR(:,:,:)< (m));

% Map containing Object without strong edges
map = map & map2;

% Estimation of the magnitude noise STD in HHH subband
Nsig = median(abs(tmp(map)))/0.6745;

% Computation of SNR on object 
fima=convn(ima,ones(3,3,3),'same');
[mu,mask]=kmeans(fima,2);
th=mean(mu);
map = find(fima>th);
SNR = mean(ima(map)) / Nsig;

% Iterative estimation of truth SNR based on Koay method
for un =1:500
    SNR2 = sqrt(epsi(SNR)*(1 + mean(ima(map))^2  / Nsig^2 )-2);
    SNR2=abs(SNR2);
    
    if abs(SNR-SNR2) < 0.000000001 
        break;
    end
    
    SNR = SNR2;
end
    

sig = sqrt((Nsig^2 / epsi(SNR)));



