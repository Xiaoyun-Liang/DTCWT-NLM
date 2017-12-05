%%%%%%%%  Main program:
%%%%%%%%  Non-local means combined with dual-tree complex wavelet transform
%%%%%%%%  for image denoising, especially arterial spin labeling (ASL)
%%%%%%%%  fMRI images
%%%%%%%%  Note: (1) The implementation of 3D dual-tree complex wavelet is mainly based on a wavelet toolbox developed by Ivan Selesnick and co-workers 
%%%%%%%%            (http://eeweb.poly.edu/iselesni/WaveletSoftware/people.html);
%%%%%%%%        (2) The implementation of onlm - optimized block-wise non-local means denoising, was released by Pierrick Coupe as a compiled mex-file 
%%%%%%%%            based on the reference: P. Coupé, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot. An Optimized Blockwise NonLocal Means 
%%%%%%%%            Denoising Filter for 3-D Magnetic Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425–441, 2008. 
%%%%%%%%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input test_data: Any images with analyze format, either single or multiple
% images
%            

% Output test_data: Denoised images are saved to a specified folder
% Two arterial spin labeling (ASL) fMRI test_datasets are included in test_test_data folder for testing
% the software
%
% Note:(1) To read and write image in analyze format, SPM and Resting-State fMRI test_data Analysis Toolkit (REST)should be
%         installed beforehand
%      (2) All related wavelet code is included in waveletpackage
%     
%
% Reference: The Matlab code is mainly based on the following paper:
%            Xiaoyun Liang, Alan Connelly, Fernando Calamante. Voxel-wise
%            functional connectomics using arterial spin labeling fMRI: the
%            role of denoising. Brain Connectivity (in press).
%
% Copyright 2015 Brain Research Institute, Melbourne, Australia
% Written by Xiaoyun Liang (Email: x.liang@brain.org.au)
% This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied  
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Add wavelet code rquiired for current implementation
clear all;
% addpath waveletpackage
% Add code for onlm and noise estimation
% addpath othercode


 % intialize the images: x=64, y=51,z=20, 59 volumes collected
 I_orig=zeros(64,51,20,59); % original images, 
 I_SP=zeros(64,51,20,59);   % denoised with samll patch, preserving features
 I_BP=zeros(64,51,20,59);   % denoised with big patch, removing noise
 I_denoised=zeros(64,51,20,59);  % denoised image by using DT-CWT-NLM
 

% Read time series
for i=1:1
    if i<10
       filename=sprintf('test_data/ASL_long_PLD/rgrase_im_single00%d.hdr',i);  % directory of input test_data
       filename1=sprintf('test_data/ASL_long_PLD/denoising_single00%d',i);    % directory of output test_data
    else
        if i<100
           filename=sprintf('test_data/ASL_long_PLD/rgrase_im_single0%d.hdr',i);
           filename1=sprintf('test_data/ASL_long_PLD/denoising_single0%d',i);
        else
           filename=sprintf('test_data/ASL_long_PLD/rgrase_im_single%d.hdr',i); 
           filename1=sprintf('test_data/ASL_long_PLD/denoising_single%d',i);
        end
    end
    
    [I_orig(:,:,:,i),header]=rest_ReadNiftiImage(filename);


    ref=I_orig(:,:,:,i);
    I_orig1=squeeze(ref);

    %estimation of noise level
    estimated_noise=RicianSTD(I_orig1);
    
    s=size(I_orig);


    % adjust the noise level to balance the effect of denoising and smoothing  
    level=estimated_noise/4; 


    % Make sure the intensity is non-negative
    L = min(ref(:));
    I_orig(:,:,:,i) = I_orig(:,:,:,i) + abs(L);

    % parameters for block-wise NLM
    M=3;    %  controlling size of search volume 
    alpha=1;  % controlling block size


    % Obtaining images with preserved featrues
    
    I_SP(:,:,:,i)=onlm(I_orig(:,:,:,i),M,alpha,level); % onlm - optimized block-wise non-local means denoising, compiled mex-file
    I_SP(:,:,:,i)=I_SP(:,:,:,i) - abs(L);


    % Obtaining images with noise components removed

    I_BP(:,:,:,i)=onlm(I_orig(:,:,:,i),M,alpha+1,level);
    I_BP(:,:,:,i)=I_BP(:,:,:,i) - abs(L);



    % Further denoising with DT-CWT
    I_denoised(:,:,:,i) = DTCWT(I_orig(:,:,:,i),I_SP(:,:,:,i),I_BP(:,:,:,i));


    % Write the denoised images to a specified folder
    rest_WriteNiftiImage(I_denoised(:,:,:,i),header,filename1);




end


