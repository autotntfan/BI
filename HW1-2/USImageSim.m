% Introduction to Biomedical Imaging   
%	MATLAB sample codes for ultrasound image simulation
%	ultrasound image = PSF convolve with scatterer distribution
%	The convolution is done in frequency domain, 
%	i.e., FT{ultrasound image} = FT{PSF} times FT{scatterer distribution}
% 	
%	The image will be with a cyst in the middle. 
%	The imaging field of view is 25.6 mm X 25.6 mm.
%
% 	user-specified parameters are:
% 		cystdB: cyst contrast relative to background in dB
%		cystR:  cyst radius in pixels
% 		displayDR: display dynamic range in dB
%		masterG: master gain
%
%	parameters for verifying the properties of speckle noises
%	EnvelopeData: envelope data
%	dBData: data with normalization and then log conversion, in dB
%
%													Revised by Meng-Lin Li, Ph. D.
%													Dept. of Electrical Engineering, National Tsing Hua University
%													11/01/2018
%													Modified by Meng-Lin Li, 10/24/2019
%													Modified by Meng-Lin Li, 10/29/2020
%													Modified by Meng-Lin Li, 03/29/2022
%

% ---------- user specified parameters ---------------
% !!! ALL PARAMETERS YOU NEED TO MODIFY ARE HERE !!!!
cystdB=+20;
cystR=50;
displayDR=60; % Dynamic range for image dispaly
masterG=0; % Gain for image display
% ----------------------------------------------------

% define raw image
xpixel=256;
ypixel=256;
pointx1=30;
pointy1=30;
pointdB=0;

cystL=10^(cystdB/20);
pointL=10^(pointdB/20);
xpixel2=(xpixel+1)/2;
ypixel2=(ypixel+1)/2;
xyaxis=([1:xpixel]'-xpixel2)*ones(1,ypixel)+sqrt(-1)*ones(xpixel,1)*([1:ypixel]-ypixel2);
mag=ones(xpixel,ypixel);
mag(find(abs(xyaxis)<cystR))=cystL*ones(size(find(abs(xyaxis)<cystR)));

randn('state',1);
rawr=randn(xpixel,ypixel);
randn('state',10);
rawi=randn(xpixel,ypixel);
raw=rawr+sqrt(-1)*rawi;
raw=raw.*mag;	% !!! scatterer distribution
raw(pointx1,pointy1)=pointL*mean(mean(abs(rawr)));


% !!! define point spread function
f0=3;
bw=1;
N=64;
d=0.25;
sint=[-ypixel/2:ypixel/2-1]/ypixel*2+eps;
faxis=[max(0.5,f0-2*bw):0.1:f0+2*bw];
tspec=exp(-pi*(faxis-f0)./bw.^2);
pa=zeros(size(sint));
for i=1:length(faxis)
   lambda=1.54/faxis(i);
	pa=pa+tspec(i)*sin(pi/lambda*N*d*sint)./sin(pi/lambda*d*sint);
end
pa=pa.^2;

psfx=fftshift(exp(-pi*(([1:xpixel]'-xpixel2)*bw/5).^2));
%psfx=(exp(-pi*(([1:xpixel]'-xpixel2)*bw/5).^2));
psfx=psfx/max(abs(psfx));
psf1=(psfx*fftshift(pa)).';
%psf1=(psfx*(pa)).';

% obtain the image
imagedata1=(ifft2(fft2(raw).*fft2(psf1))).'; % modified by Meng-Lin Li, 03/29/2022 => Now (row, column) = (z, x) (in fact, = (R, sin(theta)), sector scan format)
imagelog1=20*log10(abs(imagedata1));
imagelog1=imagelog1-max(max(imagelog1));
imageout=round((imagelog1+displayDR)*255/displayDR);
imageout(find(imageout==255))=255*ones(size(find(imageout==255)));
imageout(find(imageout<0))=zeros(size(find(imageout<0)));

% Original scatterer distribution
raw=20*log10(abs(raw));
raw=raw-max(max(raw));
raw=round((raw+displayDR)*255/displayDR);
raw(find(raw==255))=255*ones(size(find(raw==255)));
raw(find(raw<0))=zeros(size(find(raw<0)));

% display
figure
hold off
image(imagelog1+masterG+displayDR);
colormap(gray(displayDR));
colorbar
axis image

% ---------- homework related parameters ---------------
% !!! ALL PARAMETERS YOU NEED TO VERIFY THE PROPERTIES OF SPECKLE NOISES !!!!
EnvelopeData = abs(imagedata1); % Envelope Data, amplitude data; intensity = (amplitude)^2
dBData = imagelog1; % Data with normalization and then log conversion, in dB

