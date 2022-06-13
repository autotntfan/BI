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

% image parameters
nx = 128; ny = 128; %number of pixel in x, y axis
dx = 4;		                 % pixel size, in mm, 4 mm / pixel
x = dx * ([1:nx]'-(nx+1)/2);
y = -dx * ([1:ny]'-(ny+1)/2);

load ImData;
figure(1)
imagesc(ImData);   
colormap(gray)
axis image
title('Test phantom')

% ----- Create sinogram -----
% geometry parameters
nr = 128;	% # of radial samples or # of rays
dr = 4;		% ray spacing, in mm
na = ceil(pi/2*nr);		% # of views (# of angles)
r = dr * ([1:nr]'-(nr+1)/2);	% radial sample positions
angle = [0:(na-1)]'/na * pi;	% angular sample positions

% compute sinogram  
% !!! just for your reference, but not limited to it
sinogram = zeros(nr,na);
for ia = 1:na
    disp(sprintf('angle %g of %g', ia, na))
    sinogram(:,ia) = sum(imrotate(?, ?, 'bilinear', 'crop').');
end
figure(2)
imagesc(angle, r, sinogram) % !!! Please put the correct axis labels on the sinogram
title('Sinogram of the test phantom')
colorbar

%%%%%%%%%%%%%%%%%%%%%%
% QUESTION (b) -  Backprojection (Produce Laminogram), you may also refer MATLAB iradon()
%%%%%%%%%%%%%%%%%%%%%%
lamin = zeros(nr,nr); % laminogram, !!! should be zeros(nx, ny), but for simplicity, let nr = nx = ny so zeros(nx, ny) = zeros(nr, nr);
for ia = 1:na

    disp(sprintf('angle %g of %g', ia, na))
    
    tmp = imrotate(sinogram(:,ia)*ones(1,nr), ?, 'bilinear', 'crop');
    lamin = lamin+tmp;
end

figure(3)
imagesc(lamin); 
colormap(gray)
axis image
title('Laminogram')

%%%%%%%%%%%%%
% QUESTION (c) - Filter the sinogram
%%%%%%%%%%%%%
% ------ filter design
order = 256;    % zero padding, 128 to 256
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
%filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));

% --- 3) case 'cosine'
%filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));

% --- 4) case 'hamming'  
filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));

% ----5) case 'hann'
%filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;

filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the filter

% Filter Sinogram in Spatial frequency domain 
sinogramfilt = zeros(nr,na);
sinogramfilt_fft = fft(sinogram,256);  % 256 pt FFT % zero padding, 128 to 256, i.e., perform interpolation in spatial frequency domain
for ia = 1:na
   sinogramfilt_fft(:,ia) =  ?  % frequency domain filtering
end    
sinogramfilt = real(ifft(?))  % real part, "ifft()", ==> imag part should be close to zero
sinogramfilt(129:end,:) = []; % remove zero padding, from 256 pts back to 128 data pts

figure(4)
imagesc(angle, r,sinogramfilt); % !!! Please put the correct axis labels on the sinogram
colormap(gray)
xlabel('angle')
ylabel('r')
colorbar

%%%%%%%%%%%%%%%
% QUESTION (d) - Backproject the filtered sinogram
%%%%%%%%%%%%%%

ReconstructedImage = zeros(nr,nr);
for ia = 1:na
    disp(sprintf('angle %g of %g', ia, na))    
    tmp = imrotate(sinogramfilt(:,ia)*ones(1,nr), ?, 'bilinear', 'crop');
    ReconstructedImage = ReconstructedImage+tmp; 
end

% Display the reconstructed image with negative values being set to zero (reconstructed attenuation map should be >= 0)
figure(5)
imagesc(?,?,max(ReconstructedImage,0))  % please put correct axis labels on it
axis image
title()
xlabel()

%%%%%%%%%%%%%%%
% QUESTION (e) - backprojection vs. filtered backprojection
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%
% QUESTION (f) - backprojection PSF vs. filtered backprojection PSF
%%%%%%%%%%%%%%
% (1) Replace ImData with a point (i.e., only a certain pixel has non-zero attenuation coefficient.
% (2) Get the sinogram of the new ImData
% (3) Backprojection and filtered backprojection to reconstruct the PSFs

%%%%%%%%%%%%%%%
% QUESTION (g) - Effect of "View" and "Ray" on the reconstructed images and PSFs
%%%%%%%%%%%%%%



