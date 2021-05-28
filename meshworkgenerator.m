%set parameters for simulation
imageMode = 'SRRF'; %set image mode from SRRF, 3D SIM, TIRF
numFil = 25; %number of mother filaments
L =10000; %filament length in nm
%PSFs are ESTIMATES and relative to pixel size, not nm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch imageMode
    case 'SRRF'
        sizeI = 10010;
        % Size of field of view/ROI to simulate, in nm
        bin = 776;
        % number of pixels relative to sizeI
        psf = 8.02;
        % Estimated point spread function in pixels
    case '3D SIM'
        sizeI = 10020;
        bin = 310;
        psf = 3.2;
    case 'TIRF'
        sizeI = 10020;
        bin = 155;
        psf = 1.6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate blank image sized as set above
binaryImage = zeros(sizeI);

figure;
hold on;

%Generate first filaments

alpha = randi([0 360],1,1);
%randomly generated filament angle 

xS = randi([0 10000],1,1); 
yS = randi([0 10000],1,1);
%generate random seed point

x1 = xS +(L*cos(alpha));
y1 = yS +(L*sin(alpha));
%filament start point

x2 = xS -(L*cos(alpha));
y2 = yS -(L*sin(alpha));
%filament end point

hIm = imshow(binaryImage, []);
hLine = imline(gca,[x1 x2],[y1 y2]);
BIline = createMask(hLine, hIm);

%plot initial filament

    close all; 
    
LD = L; %randi([1000 3000]); %daughter filament length can be adjusted independent of mother fils
r = (abs(x1-x2)).*rand(1,1);
if x1<x2
xD = x1 + abs(r);
if y1<y2
yD = y1 + abs(tan(alpha)*r);
else 
    yD = y1 - abs(tan(alpha)*(r));
end
    
else
    xD = x2 + abs(r);
    if y1<y2
yD = y2 - abs(tan(alpha)*r);
else 
    yD = y2 + abs(tan(alpha)*(r));
end
end
%select random point on mother filament 

pos = randi([1 2]);
switch pos
    case 1
    ang = alpha + 70;
    case 2
    ang = alpha - 70;
          
end 

%chose +ve or -ve branch angle for daughter

xD2 = xD + (LD*cos(ang));
yD2 = yD + (LD*sin(ang));

% daughter filament end point

hImD = imshow(BIline, []);
hLineD = imline(gca,[xD xD2], [yD yD2]); 
BIlineD = createMask(hLineD, hImD);
BinaryMesh = imoverlay(BIline,BIlineD, 'white');

binaryImage = BinaryMesh;

%plot daughter filament

%then loops same operations for numFil set above
for i = 1:numFil
numDau = 1;  

hold on;

alpha = randi([0 360],1,1);

xS = randi([0 10000],1,1);%random seed 
yS = randi([0 10000],1,1);

x1 = xS +(L*cos(alpha));
y1 = yS +(L*sin(alpha));

x2 = xS -(L*cos(alpha));
y2 = yS -(L*sin(alpha));

hIm = imshow(binaryImage, []);
hLine = imline(gca,[x1 x2],[y1 y2]);
BIlinex = createMask(hLine, hIm);
BIline = imoverlay(binaryImage, BIlinex, 'white');

for ii = 1:numDau
    close all; 
LD = L; %randi([1000 3000]); %daughter filament length
r = (abs(x1-x2)).*rand(1,1);
if x1<x2
xD = x1 + abs(r);
if y1<y2
yD = y1 + abs(tan(alpha)*r);
else 
    yD = y1 - abs(tan(alpha)*(r));
end
    
else
    xD = x2 + abs(r);
    if y1<y2
yD = y2 - abs(tan(alpha)*r);
else 
    yD = y2 + abs(tan(alpha)*(r));
end
end
pos = randi([1 2]);
switch pos
    case 1
    ang = alpha + 70;
    case 2
    ang = alpha - 70;
          
end 
xD2 = xD + (LD*cos(ang));
yD2 = yD + (LD*sin(ang));

hImD = imshow(BIline, []);
hLineD = imline(gca,[xD xD2], [yD yD2]); 
BIlineD = createMask(hLineD, hImD);
BinaryMesh = imoverlay(BIline,BIlineD, 'white');

binaryImage = BinaryMesh;
end

end

SE = strel('square', 7); 
MeshDilate = imdilate(BinaryMesh, SE);
% apply dilation at pix per nm scale

BinnedMesh = imresize(MeshDilate, [bin bin]); 
% bin image to relevant pixel size

BinaryBinMesh = imbinarize(BinnedMesh);
% binarise image
BinaryBinMesh = squeeze(BinaryBinMesh(:,:,3));
BinaryBinMesh = im2uint8(BinaryBinMesh);
NoiseMesh = imnoise(BinaryBinMesh,'gaussian', 0.1);
% add gaussian noise
ConvMesh = imgaussfilt(NoiseMesh, psf);
% convolve with estimated PSF
SimulatedImage = imnoise(ConvMesh,'poisson');
% add poisson noise
H = fspecial('average');
SimulatedImage = imfilter(SimulatedImageN, H);
%smooth


imshow(BinnedMesh)
set(gcf,'Position', [0 0 bin bin]);
export_fig groundtruthMesh.tif;

imshow(SimulatedImage)
set(gcf,'Position', [0 0 bin bin]);
export_fig convMesh.tif;
