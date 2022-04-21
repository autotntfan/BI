% EE441000 林祐安 110061622 HW1 part1 04/09/2022
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
dz = z_axis(2)-z_axis(1);
x_axis = (0:1:(n-1))*dx - (n-1)*dx/2;    % x axis
envelope = abs(hilbert(rf_data));   % envelope detection
envelope_dB = 20*log10(envelope/max(max(envelope))+eps);    % log conversion with respect to the maximum value 
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
% PSF_x = max(PSF); % lateral beam profile. projection along the depth. 
% Give it a thought: RF PSF or envelope-detected PSF, which one should you use? 
% Remember to elaborate in your report.
% image_x = conv(scatterer_x, PSF_x); % vary the distance, and then try the eye examination
Distance = 0.45; % in mm
Distance = Distance/dz;
scatterer = zeros(m,n);
scatterer(round((m-Distance)/2),round(n/2)) = 1;
scatterer(round((m+Distance)/2),round(n/2)) = 1;
axial_img = zeros(1,length(pt_location));
figure
for ii = 1:length(pt_location)
    RF_PSF = rf_data(floor(1+(ii-1)*m/length(pt_location)):floor(ii*m/length(pt_location)),:);
    img = conv2(scatterer,RF_PSF,'same');
    envelope_img = abs(hilbert(img));
    envelope_img_dB = 20*log10(envelope_img/max(max(envelope_img))+eps);
    subplot(1,5,ii)
    axial_project = max(envelope_img,[],2);
    axial_img(ii) = sum(axial_project > (max(axial_project)-6)) * dz;
    image(envelope_img_dB + DR)
    colormap(gray(DR))
    if ii == 5
        Far_PSF = RF_PSF;
    elseif ii == 3
        Focal_PSF = RF_PSF;
    end
end

Distance = 0.75; % in mm
Distance = Distance/dx; 
scatterer = zeros(m,n);
scatterer(round(m/2),round((n-Distance)/2)) = 1;
scatterer(round(m/2),round((n+Distance)/2)) = 1;
lateral_img = zeros(1,length(pt_location));
figure
for ii = 1:length(pt_location)
    RF_PSF = rf_data(floor(1+(ii-1)*m/length(pt_location)):floor(ii*m/length(pt_location)),:);
    img = conv2(scatterer,RF_PSF,'same');
    envelope_img = abs(hilbert(img));
    envelope_img_dB = 20*log10(envelope_img/max(max(envelope_img))+eps);
    subplot(1,5,ii)
    lateral_project = max(envelope_img,[],1);
    lateral_img(ii) = sum(lateral_project > (max(lateral_project)-6)) * dx;
    image(envelope_img_dB + DR)
    colormap(gray(DR))
end

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


lateral_res6dB = zeros(1,length(pt_location));
axial_res6dB = zeros(1,length(pt_location));
lateral_res20dB = zeros(1,length(pt_location));
axial_res20dB = zeros(1,length(pt_location));
for ii = 1:length(pt_location)
    patch = envelope_dB(floor(1+(ii-1)*m/length(pt_location)):floor(ii*m/length(pt_location)),:);
    lateral_project = max(patch,[],1);
    lateral_res6dB(ii) = sum(lateral_project > (max(lateral_project)-6)) * dx;
    lateral_res20dB(ii) = sum(lateral_project > (max(lateral_project)-20)) * dx;
    axial_project = max(patch,[],2);
    axial_res6dB(ii) = sum(axial_project > (max(axial_project)-6)) * dz;
    axial_res20dB(ii) = sum(axial_project > (max(axial_project)-20)) * dz;
end


% (f)
% Find the aperture size of the transducer based on the theoretic lateral resolution "at the focal point"
ALINE = abs(fftshift(fft(rf_data(:,round(n/2)))));
freq_axis = linspace(-fs/2,fs/2,length(ALINE));
[~, ind] = max(ALINE); % center freq. of the transducer, in MHz, try "fft" the provided RF data to find out the center frequency. Give it a thought: Which A-line and which point-target data should you use? Again, remember to elaborate this point in your report.
fc = abs(freq_axis(ind));
lambda = soundv/fc; % wavelength of the center frequency
depth = focal_pt;
f_number = lateral_res6dB(3)/lambda;
aper_size = focal_pt/f_number;  % aperture size, in mm
figure
plot(freq_axis,ALINE)
xlabel('freq. (MHz)')


% (g)

% (h)

% (i)
scale = 40; % in mm
N = 1e5; % number of scatterers
Nz = round(scale/dz);
Nx = round(scale/dx);
scatterer_pos_x = rand(N, 1)*scale; % scale: used to scale the value to within 4 cm, x position, position unit or lateral beam spacing (i.e., lateral spatial sampling interval) has to be the same as that of PSF
scatterer_pos_z = rand(N, 1)*scale; % scale: used to scale the value to within 4 cm, z position, position unit or axial sampling interval has to be the same as that of PSF
% or
% scatter_pos = rand(N)*scale;
scatterer_dist = zeros(Nz,Nx); % spatial distribution of the scatterers, 
                               % Nz: grid points along the z direction, determined by the "4 cm" field of view in z and the spatial sampling interval in the z direction
							   % Nx: grid points along the x direction, determined by the "4 cm" field of view in x and the spatial sampling interval in the x direction

scatterer_pos_x = ceil(scatterer_pos_x/dx);
scatterer_pos_z = ceil(scatterer_pos_z/dz);
% locate the scatterers into the scatterer distribution matrix
for ii = 1:N
    iX = scatterer_pos_x(ii);
    iZ = scatterer_pos_z(ii);
    if ((iX*dx-20)^2+(iZ*dz-20)^2) <= 5^2
        scatterer_dist(iZ, iX) = scatterer_dist(iZ, iX) + 0.1;
    else
        scatterer_dist(iZ, iX) = scatterer_dist(iZ, iX) + 1;
    end
end

Focal_image = conv2(scatterer_dist, Focal_PSF, 'same'); % note that the sampling interval (in x and in z) of the scatterer distribution should be the same as that of RF_PSF.
                                          % RF_PSF: PSF in RF data form
										  % RF_image: RF data
Far_image = conv2(scatterer_dist, Far_PSF, 'same');
% then envelope detection with Hilbert transform, log conversion, and then determine the image dynamic range to show the image 
% then estimate the contrast
envelope_img = abs(hilbert(Focal_image));
dB_img = 20*log10(envelope_img/max(max(envelope_img))+eps);
figure
image(dB_img + DR)
colormap(gray(DR))
title('at focal zone')

envelope_img = abs(hilbert(Far_image));
dB_img = 20*log10(envelope_img/max(max(envelope_img))+eps);
figure
image(dB_img + DR)
colormap(gray(DR))
title('at 18 mm')

% calculate contrast ratio
ROI = zeros(Nz,Nx);
ROI(round(17/dz):round(23/dz),round(17/dx):round(23/dx)) = 1;
ROI_img = ROI.*dB_img;
CR_ROI = sum(ROI_img,'all')/sum(ROI,'all'); % average contrast in ROI
BKG = zeros(Nz,Nx);
BKG(round(1/dz):round(10/dz),round(1/dx):round(10/dx)) = 1;
BKG_img = BKG.*dB_img;
CR_BKG = sum(BKG_img,'all')/sum(BKG,'all'); % average contrast in background
CR = CR_ROI - CR_BKG;

% (j) Bonus
DemodData = rf_data.*repmat(exp(-1j*2*pi*fc*time_axis'),1,n); % freq. shift
f_order = 40;
fcut = 0.2;
f = [0 fcut fcut 1];
mag = [1 1 0 0];
b1 = fir2(f_order,f,mag); % design a LPF
BaseData = abs(filter(b1,1,DemodData));
envelope_dB = 20*log10(BaseData/max(max(BaseData))+eps);
figure
image(x_axis,z_axis, envelope_dB + DR)
colormap(gray(DR))
axis image