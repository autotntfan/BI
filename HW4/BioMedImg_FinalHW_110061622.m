% Introduction to Biomedical Imaging,   Spring 2022
%   Template of HW4
%	See also iradon() in MATLAB
%					 				Edited by Meng-Lin Li, Ph. D., 12/27/2016
%									Dept. of Electrical Engineering, National Tsing Hua University
%
%									Revised by Meng-Lin Li, Ph. D., 01/04/2018
%									Revised by Meng-Lin Li, Ph. D., 12/24/2020
%									Revised by Meng-Lin Li, Ph. D., 05/24/2022
%


%%%%%%%%%%%%%%
% Question (a)
%%%%%%%%%%%%%%
clear 
close all

load ImData;
% ImData = ImData(1:2:end,1:2:end);
% image parameters
[nx, ny] = size(ImData); %number of pixel in x, y axis
dx = 4;		                 % pixel size, in mm, 4 mm / pixel
x = dx * ([1:nx]'-(nx+1)/2);
y = -dx * ([1:ny]'-(ny+1)/2);


figure(1)
imagesc(x,-y,ImData);   
colormap(gray)
axis image
title('Test phantom')

% ----- Create sinogram -----
% geometry parameters
nr = nx;	% # of radial samples or # of rays
dr = 4;		% ray spacing, in mm
na = ceil(pi/2*nr);		% # of views (# of angles)
r = dr * ([1:nr]'-(nr+1)/2);	% radial sample positions
angle = [0:(na-1)]'/na * pi;	% angular sample positions

% compute sinogram  
% !!! just for your reference, but not limited to it
sinogram = zeros(nr,na);
for ia = 1:na
    sinogram(:,ia) = sum(imrotate(ImData, rad2deg(angle(ia)), 'bilinear', 'crop')');
end
figure(2)
imagesc(angle, r, sinogram) % !!! Please put the correct axis labels on the sinogram
title('Sinogram of the test phantom')
ylabel('radius')
xlabel('angular')
colorbar

%%%%%%%%%%%%%%%%%%%%%%
% QUESTION (b) -  Backprojection (Produce Laminogram), you may also refer MATLAB iradon()
%%%%%%%%%%%%%%%%%%%%%%
lamin = zeros(nr,nr); % laminogram, !!! should be zeros(nx, ny), but for simplicity, let nr = nx = ny so zeros(nx, ny) = zeros(nr, nr);
for ia = 1:na 
    tmp = imrotate(sinogram(:,ia)*ones(1,nr), -rad2deg(angle(ia)), 'bilinear', 'crop');
    lamin = lamin+tmp;
end
figure(3)
imagesc(x,-y,lamin); 
colormap(gray)
axis image
title('Laminogram')

%%%%%%%%%%%%%
% QUESTION (c) - Filter the sinogram
%%%%%%%%%%%%%
% ------ filter design
order = 2*nx;    % zero padding, 128 to 256
d = 1;
% First create a ramp filter - go up to the next highest
% power of 2.
filt = 2*( 0:(order/2) )./order;    % filter
w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist 

% !!! uncomment the filter you choose
% --- 1) case 'ram-lak'
% Do nothing

% --- 2) case 'shepp-logan'
% be careful not to divide by 0:
% filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));

% --- 3) case 'cosine'
% filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));

% --- 4) case 'hamming'  
% filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));

% ----5) case 'hann'
% filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;

filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the filter

% Filter Sinogram in Spatial frequency domain 
sinogramfilt = zeros(nr,na);
sinogramfilt_fft = fft(sinogram,2*nx);  % 256 pt FFT % zero padding, 128 to 256, i.e., perform interpolation in spatial frequency domain
for ia = 1:na
   sinogramfilt_fft(:,ia) =  sinogramfilt_fft(:,ia).*filt;  % frequency domain filtering
end    
sinogramfilt = real(ifft(sinogramfilt_fft));  % real part, "ifft()", ==> imag part should be close to zero
sinogramfilt(nx+1:end,:) = []; % remove zero padding, from 256 pts back to 128 data pts

figure(4)
imagesc(angle, r,sinogramfilt); % !!! Please put the correct axis labels on the sinogram
xlabel('angle')
ylabel('r')
colorbar

%%%%%%%%%%%%%%%
% QUESTION (d) - Backproject the filtered sinogram
%%%%%%%%%%%%%%

ReconstructedImage = zeros(nr,nr);
for ia = 1:na
    tmp = imrotate(sinogramfilt(:,ia)*ones(1,nr), -rad2deg(angle(ia)), 'bilinear', 'crop');
    ReconstructedImage = ReconstructedImage+tmp; 
end
ReconstructedImage = max(ReconstructedImage,0);
% Display the reconstructed image with negative values being set to zero (reconstructed attenuation map should be >= 0)
figure(5)
imagesc(ReconstructedImage)  % please put correct axis labels on it
colormap(gray)
axis image
title('Reconstruction Image')
xlabel('position')
hold on
lines = 25;
plot([1 nx],[lines lines])
%%%%%%%%%%%%%%%
% QUESTION (e) - backprojection vs. filtered backprojection
%%%%%%%%%%%%%%
Aline_original = 20*log10(lamin/max(max(lamin))+eps);
Aline_original = Aline_original(lines,:);
Aline_filt = 20*log10(ReconstructedImage/max(max(ReconstructedImage))+eps);
Aline_filt = Aline_filt(lines,:);
figure
plot(Aline_original,'black')
hold on
plot(Aline_filt,'red')
legend('orginal','filtered')
FWHM_original = sum(Aline_original>max(Aline_original)-6)*dx;
FWHM_filt = sum(Aline_filt>max(Aline_filt)-6)*dx;
%%%%%%%%%%%%%%%
% QUESTION (f) - backprojection PSF vs. filtered backprojection PSF
%%%%%%%%%%%%%%
% (1) Replace ImData with a point (i.e., only a certain pixel has non-zero attenuation coefficient.
% (2) Get the sinogram of the new ImData
% (3) Backprojection and filtered backprojection to reconstruct the PSFs
ImData = zeros(size(ImData));
ImData(nx/2,nx/2) = 1;
sinogram = zeros(nr,na);
for ia = 1:na
    sinogram(:,ia) = sum(imrotate(ImData, rad2deg(angle(ia)), 'bilinear', 'crop')');
end
figure
imagesc(angle, r, sinogram) % !!! Please put the correct axis labels on the sinogram
title('Sinogram of the test phantom')
ylabel('radius')
xlabel('angular')
colorbar

lamin = zeros(nr,nr); % laminogram, !!! should be zeros(nx, ny), but for simplicity, let nr = nx = ny so zeros(nx, ny) = zeros(nr, nr);
for ia = 1:na 
    tmp = imrotate(sinogram(:,ia)*ones(1,nr), -rad2deg(angle(ia)), 'bilinear', 'crop');
    lamin = lamin+tmp;
end
figure
imagesc(lamin); 
colormap(gray)
axis image
title('Laminogram')

% Filter Sinogram in Spatial frequency domain 
sinogramfilt = zeros(nr,na);
sinogramfilt_fft = fft(sinogram,2*nx);  % 256 pt FFT % zero padding, 128 to 256, i.e., perform interpolation in spatial frequency domain
for ia = 1:na
   sinogramfilt_fft(:,ia) =  sinogramfilt_fft(:,ia).*filt;  % frequency domain filtering
end    
sinogramfilt = real(ifft(sinogramfilt_fft));  % real part, "ifft()", ==> imag part should be close to zero
sinogramfilt(nx+1:end,:) = []; % remove zero padding, from 256 pts back to 128 data pts

figure
imagesc(angle, r,sinogramfilt); % !!! Please put the correct axis labels on the sinogram
xlabel('angle')
ylabel('r')
colorbar



ReconstructedImage = zeros(nr,nr);
for ia = 1:na
    tmp = imrotate(sinogramfilt(:,ia)*ones(1,nr), -rad2deg(angle(ia)), 'bilinear', 'crop');
    ReconstructedImage = ReconstructedImage+tmp; 
end
ReconstructedImage = max(ReconstructedImage,0);
% Display the reconstructed image with negative values being set to zero (reconstructed attenuation map should be >= 0)
figure
imagesc(ReconstructedImage)  % please put correct axis labels on it
colormap(gray)
axis image
title('Reconstruction Image')
xlabel('position')
lines = nx/2;
% lamin = interp2(lamin,2);
% ReconstructedImage = interp2(ReconstructedImage,2);
Aline_original = 20*log10(lamin/max(max(lamin))+eps);
Aline_original = Aline_original(lines,:);
Aline_filt = 20*log10(ReconstructedImage/max(max(ReconstructedImage))+eps);
Aline_filt = Aline_filt(lines,:);
figure
plot(Aline_original,'black')
hold on
plot(Aline_filt,'red')
legend('orginal','filtered')
FWHM_original = sum(Aline_original>max(Aline_original)-6)*dx;
FWHM_filt = sum(Aline_filt>max(Aline_filt)-6)*dx;


%%%%%%%%%%%%%%%


% QUESTION (g) - Effect of "View" and "Ray" on the reconstructed images and PSFs
%%%%%%%%%%%%%%



