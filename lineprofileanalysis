I = imread('image to analyse');
minheight = 50; %set min peak height - sample relative
avpkdist = ["Average distance btw peaks"; ];
numpks = ["Number of peaks"; ];
linepos = ["Line angle"; ];
allwidth = ["Mean peak widths"; ];


for ii = 1:100

xS = size(I,2)/2;
yS = size(I,2)/2;
L = yS;

alpha = randi([0 360],1,1);
ang = randi(2);

x1 = xS +(L*cos(alpha));
y1 = yS +(L*sin(alpha));

  
x2 = xS -(L*cos(alpha));
y2 = yS -(L*sin(alpha));  

    
x = [x1 x2];
y = [y1 y2];

linepos = [linepos; alpha];
%Generate line across image to be analysed at random angle

c = improfile(I,x,y);
cMin = min(c);
cNorm = c - cMin;
% plot intensity profile across line and normalise

[pks,locs,widths,proms] = findpeaks(cNorm,'MinPeakHeight', minheight, 'WidthReference','halfheight'); 
%identify peaks
pkdist = mean(diff(locs));
meanwidth = mean(widths);
number = size(pks, 1);
%calculate outputs
avpkdist = [avpkdist; pkdist];
numpks = [numpks; number];
allwidth = [allwidth; meanwidth];

end

outputs = [linepos, numpks, avpkdist, allwidth];
