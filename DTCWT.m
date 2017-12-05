% Reference: The Matlab code is mainly based on the following paper:
%            Xiaoyun Liang, Alan Connelly, Fernando Calamante. Voxel-wise
%            functional connectomics using arterial spin labeling fMRI: the
%            role of denoising. Brain Connectivity (in press).

function I3 = DTCWT(I0,I1,I2)

% Denoising by using Dual-tree complex wavelet (DT-CWT)

% Input:
%     I0- original noisy 3D images
%     I1-denoising with small patch-feature preserved
%     I2-denoising with big patch-noise removed
% Output
%     I3-denoised image by using DT-CWT
% 

s = size(I1);

L1 = 2^(ceil(log2(s(1))));
L2 = 2^(ceil(log2(s(2))));
L3 = 2^(ceil(log2(s(3))));

x1 = zeros(L1,L2,L3);
x2=x1;
x3=x1;
x1(1:s(1),1:s(2),1:s(3)) = I1(:,:,:);
x2(1:s(1),1:s(2),1:s(3)) = I2(:,:,:); 
x3(1:s(1),1:s(2),1:s(3)) = I0(:,:,:); 

  J = 3;    %J=3
  [Faf, Fsf] = FSfarras;
  [af, sf] = dualfilt1;
  w1 = cplxdual3D(x1, J, Faf, af);
  w2 = cplxdual3D(x2, J, Faf, af);
  w3 = cplxdual3D(x3, J, Faf, af);
  I=sqrt(-1);
  

% wavelet thresholding using Bayes Shrink
for j=1:J
    for m=1:2
        for n=1:2
            for p=1:2
                for i=1:7
                    
  
    
                    W = w3{J}{m}{n}{p}{i};
                    
                    
                    Nv=(median(abs(W(:)))/0.6745)^2;  % estimation of the noise variance Nv
                    % standard deviation of the subband of noisy image
                    SIGMAy = std(W(:));
                    % standard deviation of the original image 
                    SIGMAx = SIGMAy^2 - Nv;
                    
                    
                    % The bayes shrink thresholding employs soft
                    % thresholding with adaptive data-driven near optimal
                    % threshold
                    if (SIGMAx < 0)
                        Tbs = max(abs(w3{j}{m}{n}{p}{i}(:)));
                    else
                        Tbs = Nv / sqrt(SIGMAx);

                    end;
                   

                    w3{j}{m}{n}{p}{i} = abs(w3{j}{m}{n}{p}{i});

                    dist = w3{j}{m}{n}{p}{i} - Tbs;
                    dist = exp(-0.01*dist);
                    dist = 1./(1+dist);

                    w3{j}{m}{n}{p}{i} =  dist.*w1{j}{m}{n}{p}{i} + (1-dist).*w2{j}{m}{n}{p}{i};

                end
            end
        end
    end

end

 
w3{J+1} = w1{J+1};

% Transform back to image space
I3 =  icplxdual3D(w3,J,Fsf,sf);

I3 = I3(1:s(1),1:s(2),1:s(3));
