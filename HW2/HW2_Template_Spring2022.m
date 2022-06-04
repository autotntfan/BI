% Introduction to Biomedical Imaging,   Spring 2022
%   Template of HW2
%					Edited by Meng-Lin Li, Ph. D., 11/28/2018
%					Dept. of Electrical Engineering, National Tsing Hua University
%					Modified by Meng-Lin Li, Ph. D., 11/20/2019
%					Modified by Meng-Lin Li, Ph. D., 11/25/2020
%					Modified by Meng-Lin Li, Ph. D., 05/05/2022
%

% ----------- Quiz 1 ------------
% --- read the provided data
% read the real component
fid = fopen('hw2_r.dat','rb');
data_r = fread(fid,[256 256],'int32').'; % why .' ?, check the reconstructed images with and without .', you will find it out
fclose(fid);
% read the imaginary component
fid = fopen('hw2_i.dat','rb');
data_i = fread(fid,[256 256],'int32').';
fclose(fid);

KspaceData = data_r+sqrt(-1)*data_i; % k-space data !!!


% --- (a) ---
% Show the real and imaginary components, respectively (by Matlab instruction: imagesc()), 
% and check the data range (by Matlab instructions: colorbar, max() and min(), or mean2()... etc.).
% -- real part
figure
imagesc(???);
axis image
axis off
colorbar

% -- imaginary part
figure
imagesc(???);
axis image
axis off
colorbar

% Do you find out that the data suffer a constant DC offset?

% --- (b) ---
% Discuss the effects on the image when the receiver channel has a constant DC offset arising from the electronics
% Justify your findings
ImData = abs(ifftshift(ifft2(???))); % image data !!!
figure
imagesc(???)
axis image
axis off
colormap(gray)

% --- (c) ---
% remove DC offset
KspaceData_DCRemoved = ???; % k-space data with DC offset being removed
DCRemovedImData = abs(ifftshift(ifft2(???)));
figure
imagesc(???)
axis image
axis off
colormap(gray)

% --- (d) ---
% Reconstruct the image by using only the even samples along the ky direction 
NewKspaceData = KspaceData_DCRemoved(???); 		 % !!! e.g., NewKspaceData = KspaceData_DCRemoved(2:2:end,:);  
												 % Note that the above codes I provide simply show you how to use only even samples in certain k-space direction. 
                                                 % Make sure you are using even samples in ky direction.
												 % If you have troubles in removing the DC offset in the original k-space data, you can load DC removed k space data from the provided file - KspaceData_DCremoved.mat
NewImData = abs(ifftshift(ifft2(???)));
figure
imagesc(???)
axis image
axis off
colormap(gray)

% --- (e) ---
% Reconstructe the image with the middle (#129) and the first (#1) phase encoding values being set to zero, respectively
KspaceData_PhaseEncodingMissed = KspaceData_DCRemoved;
KspaceData_PhaseEncodingMissed(???) = 0; % set the first or the middle phase encoding value to zero.
ImData_PhaseEncodingMissed = abs(ifftshift(ifft2(???)));
figure
imagesc(???)
axis image
axis off
colormap(gray)

% --- (f) ---
% Add EMI to k-space data, and 2D inverse Fourier transform the data
EMIKspaceData = KKspaceData_DCRemoved;
EMIKspaceData(100,100) = ???; % ???: a large enough number, can be real, can be complex (real, complex, any diffences?)
EMIImData = abs(ifftshift(ifft2(???)));
figure
imagesc(???)
axis image
axis off
colormap(gray)

% --- (g) ---
% Use the same set of data to try implementing half Fourier imaging which saves you about half of the scan time. 
% acquire half of k-space data, remove DC offset, 
%then use the conjugate-symmetry property of Fourier transform for a real signal to build the full k-space and reconstruct the MRI image
HalfFourierKsapceData = KspaceData_DCRemoved(1:128,:); % you may try 1:128, or 1:129
ReconstructedCompleteKspaceData = ???
HalfFourierImData = abs(ifftshift(ifft2(???)));
figure
imagesc(???)
axis image
axis off
colormap(gray)

% --- (h) ---
% Based on the image you reconstruct in (c), please tell what weighted image (e.g., T1 weighted or T2 weighted) it is, and justify your answer.

% --- (i) ---
% A technician has reversed the wires to the slice selection gradient coil such that the gradient is reversed. Must you change the current waveform delivered to this gradient coil for slice selection as a result of this error? Explain 
% That is, What parameters of the current waveform, e.g., frequency and pulse duration do you need to change for the same slice selection?



% ----------- Quiz 2 ------------
% Hand writing











