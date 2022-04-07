% Introduction to Biomedical Imaging,   Fall 2020
%   Template of HW1, Part1
%                                       Edited by Meng-Lin Li, 10/18/2018
%										Modified by Meng-Lin Li, 10/15/2019
%										Modified by Meng-Lin Li, 10/19/2020
%										Modified by Meng-Lin Li, 03/17/2022
%												    Dept. of Electrical Engineering
%													National Tsing Hua University, Taiwan
clear
close all

% Note that you can use "SI" unit instead of the units I used below
time_offset = 6.48; % in usec
fs = 50;    % sampling rate, in MHz
% aper_size = ??;  % aperture size, in mm
focal_pt = 12;   % focal point, in mm
dx = 0.050;    % distance between two successive scanning positions, in mm
soundv = 1.5;  % speed of sound, in mm/us
pt_location = [6 9 12 15 18]; % points location in mm

% (a)
DR = 45;    % dynamic range, in dB
load('points_rf_data');
rf_data = points_rf_data; 
clear points_rf_data;
[m,n] = size(rf_data);
time_axis = time_offset + (0:1:(m-1))*(1/fs); % in usec
z_axis = soundv*time_axis./2;   % z axis
x_axis = (0:1:(n-1))*dx;    % x axis
envelope = abs(hilbert(rf_data));   % envelope detection
envelope_dB = 20*log(envelope/max(max(envelope))+eps);    % log conversion with respect to the maximum value 
figure
image(x_axis, z_axis, envelope_dB+DR); % or imagesc()?
colormap(gray(DR))
colorbar;
title('Point targets')
xlabel('Lateral position (mm)')
ylabel('Depth (mm)')
axis image

% (b)
% Observe the image in (a) and find the focal length of the transducer
% one of the point targets is located at the focal point

% (c)
% Based on your answer in (b), tell the speed of sound
% soundv = ??;  % speed of sound, in mm/us

% (d)
% Check the spatial resolution based on the fundamental definition of the spatial resolution and your eye examination
% for example, in 1D (you may try to do the following procedure in 2D) 
% distance = 10; % in terms of sample points. => you have to convert true distance to sample points
% scatterer_x = [ zeros(1,100) 1 zeros(1, distance) 1 zeros(1,100)] % scatterer positions in x direction
% PSF_x = max(PSF); % lateral beam profile. projection along the depth. Give it a thought: RF PSF or envelope-detected PSF, which one should you use? Remember to elaborate in your report.
% image_x = conv(scatterer_x, PSF_x); % vary the distance, and then try the eye examination


%(e)
% Check the -6ddB and -20 dB lateral and axial resolution point target by point target
% e.g., for point target 1, PSF_point1 = envelope_dB(1:190,:);
% then, find the FWHM (-6dB) and - 20dB for the normalized axial and lateral beam plot,
% respectively.

% Find the FWHM?
% Hint: quick solution by using Matlab built in function "find()"
%       but "not accurate enough"
%	for better accuracy, you need to interpolate the provided data (e.g., by "interp()" or by "resample()")
%   you may also try your own codes
% for example,
% idx = find(PSF >= 0.5); % PSF is the projected PSF after interpolation
% FWHM_PSF = (idx(end) - idx(1))*spatial sampling interval;  % in um

lateral_project = max(envelope,[],2);
figure
for ii = 1:length(pt_location)
    PSF = lateral_project(floor(1+(ii-1)*m/length(pt_location)):floor(ii*m/length(pt_location)));
    lateral_res = sum(PSF > 0.5*max(PSF)) * dz;
    subplot(1,5,ii)
    plot(project)
end
figure
plot(z_axis,lateral_project)

axial_project = max(envelope,[],1);
figure
plot(x_axis, axial_project)
% 
% % (f)
% % Find the aperture size of the transducer based on the theoretic lateral resolution "at the focal point"
% fc = ?; % center freq. of the transducer, in MHz, try "fft" the provided RF data to find out the center frequency. Give it a thought: Which A-line and which point-target data should you use? Again, remember to elaborate this point in your report.
% lambda = ?; % wavelength of the center frequency
% depth = ?;
% f_number = ?;
% aper_size = ??;  % aperture size, in mm
% 
% 
% % (g)
% 
% % (h)
% 
% % (i)
% N = ??; % number of scatterers
% scatterer_pos_x = rand(N, 1)*scale; % scale: used to scale the value to within 4 cm, x position, position unit or lateral beam spacing (i.e., lateral spatial sampling interval) has to be the same as that of PSF
% scatterer_pos_z = rand(N, 1)*scale; % scale: used to scale the value to within 4 cm, z position, position unit or axial sampling interval has to be the same as that of PSF
% % or
% % scatter_pos = rand(N)*scale;
% scatterer_dist = zeros(Nz,Nx); % spatial distribution of the scatterers, 
%                                % Nz: grid points along the z direction, determined by the "4 cm" field of view in z and the spatial sampling interval in the z direction
% 							   % Nx: grid points along the x direction, determined by the "4 cm" field of view in x and the spatial sampling interval in the x direction
% 
% % locate the scatterers into the scatterer distribution matrix
% for iX = 1:Nx,
% 	for iZ = 1:Nz,
% 		if ??? % if the grid matches the scatterer position (scatter_pos_x, scatter_pos_z)
% 			scatterer_dist(iZ, iX) = 1; % assign the same back-scattering coef.
% 		end
% 	end
% end
% 
% RF_image = conv2(scatterer_dist, RF_PSF); % note that the sampling interval (in x and in z) of the scatterer distribution should be the same as that of RF_PSF.
%                                           % RF_PSF: PSF in RF data form
% 										  % RF_image: RF data
% % then envelope detection with Hilbert transform, log conversion, and then determine the image dynamic range to show the image 
% % then estimate the contrast
% 
% 
% % (j) Bonus
% 
