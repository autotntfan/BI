% EE441000 林祐安 110061622 HW1 PartII 04/19/2022
% Introduction to Biomedical Imaging,   Spring 2022
%   Template of HW1, PartII - Expericing Ultrasound Speckle Noises by Simulations
%	
%                                       Edited by Meng-Lin Li, 11/01/2018
%										Modified by Meng-Lin Li, 10/24/2019
%										Modified by Meng-Lin Li, 10/29/2020
%										Modified by Meng-Lin Li, 03/29/2022
%										Dept. of Electrical Engineering,										
%										National Tsing Hua University

% ---------- run the provided simulation script of ultrasound B-mode images ----------
% to complete your HW, You need to modify "user-specified parameters", e.g., "cystdB" and "cystR" in USImageSim.m, 
% and you will need the variables "EnvelopeData" and "dBData" to verify the properties of speckle noises
clear
close all

USImageSim;
% (a)
% try hist()
CystRegion = EnvelopeData(round(128-cystR/2):round(128+cystR/2),round(128-cystR/2):round(128+cystR/2));
CystRegiondB = dBData(round(128-cystR/2):round(128+cystR/2),round(128-cystR/2):round(128+cystR/2));
CystIntensity = CystRegion.^2;

BKGRegion = EnvelopeData(1:round(128/4),1: round(128/4));
BKGRegiondB = dBData(1:round(128/4),1:round(128/4));
BKGIntensity = BKGRegion.^2;

figure
subplot(1,2,1)
histogram(CystRegion,'Normalization','probability')
title('magnitude histogram in cyst region')
set(gca,'FontSize',16)
subplot(1,2,2)
histogram(CystIntensity,'Normalization','probability')
title('intensity histogram in cyst region')
set(gca,'FontSize',16)

figure
subplot(1,2,1)
histogram(BKGRegion,'Normalization','probability')
title('magnitude histogram in background region')
set(gca,'FontSize',16)
subplot(1,2,2)
histogram(BKGIntensity,'Normalization','probability')
title('intensity histogram in background region')
set(gca,'FontSize',16)

% (b) 
% try mean2() and std2() for SNRI (intensity) and SNRE (amplitude, i.e., envelope data) estimation in the speckle background and the inclusion, respectively
% what's the difference between std() ( or std(std()) ) and std2(), and mean() ( or mean(mean()) ) and mean2()?
SNRE_background = mean(BKGRegion,'all')/std2(BKGRegion); % envelope data: EnvelopeData,  Inclusion: EnvelopeData(105:155,110:155), background: EnvelopeData(190:250, 50:200)
SNRE_inclusion = mean(CystRegion,'all')/std2(CystRegion);

SNRI_background = mean(BKGIntensity,'all')/std2(BKGIntensity); % Intensity = (EnvelopeData).^2;
SNRI_inclusion = mean(CystIntensity,'all')/std2(CystIntensity);

% (c)try std2(), what's the difference between std() ( or std(std()) ) and std2()?, std() got to be used along with reshape() in this case
std_inclusion_dB = std2(CystRegiondB);  % ~ close to 4.34 dB?
std_background_dB = std2(BKGRegiondB);  

% (d) change the value of the variable “cystdB” in USImageSim.m to see by your naked eyes what the minimum detectable contrast of the ultrasound imaging system is. 

% (e) change the diameter of the higher or lower scattering inclusion (i.e., change the value of the variable "cystR" in USImageSim.m) to see if you can detect better the inclusion with a larger diameter

% (f) Apply N by N moving average filter, which is a low pass filter, to log compressed data (i.e., dBData) to reduce the speckle noises
%	Can you see the contrast or restore the original contrast now when the original contrast in between the inclusion and background smaller than the minimum detectable contrast of the ultrasound imaging system?
%	If you can see the contrast now, What's the price you pay for?
N = 5; % filter/kernel size of the moving average filter, in terms of pixel number
MVF = ones(N,N)/(N*N); % N by N moving average filter
LPFData = conv2(dBData, MVF,'valid'); % low pass filtered data, Give it a thought. Why I use 'valid' instead of 'same' or 'full'?
									  % Note that, in this case, supposedly, you should barely see or not be able to see the inclusion in dBData.

% Display the low pass filtered data by imagesc()
% Can you see the inclusion now?
figure
imagesc(LPFData+displayDR)
colormap(gray)
colorbar
axis image

% Plot and compare the lateral profiles across the background and inclusion with and without low pass filtering
% Does the low pass filtering restore the original contrast?
% From the lateral profile with low pass filtering, do you know what price you pay for to restore/"see" the contrast?
BeforeLPF = dBData(128,:);
AfterLPF = [0 0 LPFData(128,:) 0 0]; % align the length 
figure
plot(BeforeLPF)
hold on
plot(AfterLPF)
legend('before','after')
ylabel('dB')






