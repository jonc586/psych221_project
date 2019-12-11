ieInit;
light_mus = 400:4:476; %provide a set of light spectrums to iterate over (this array describes mu)
light_buffer = ones(1,20);
contrast_metric = zeros(101,31);
lowpass_cutoff = linspace(380,780,101);
highpass_cutoff = linspace(400,700,31);
for b=1 %choose which light mus you wanna iterate over, can go from 1:20
%% Parameterize light and fluorescence
waveA = 380:4:780;
deltaL = waveA(2) - waveA(1);
nWaves = length(waveA);
% Read in ISETCam stock lighting options
fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,waveA); % 14 lights
nChannels = size(illuminant,2);
illuminantPhotons = Energy2Quanta(waveA,illuminant);
whichLight = 2;    % select a light you wanna use
sceneIlluminant = illuminant(:,whichLight)';
%plot(linspace(380,780,101),sceneIlluminant);
%Override the ISETCam stock lighting options if you wanna use a light we
%specify with a gaussian. Comment this part out if you're not interested in
%doing this and just wanna use stock options
x = linspace(380,780,101);
sceneIlluminant = normpdf(x,light_mus(b),35);
sceneIlluminant = sceneIlluminant./max(sceneIlluminant);
sceneIlluminant = sceneIlluminant.* 0.07; 
nPixels = 64; %number of pixels you want in the scene
%Create a square tile where we can apply fluorophores
scene = sceneCreate('uniform equal energy',nPixels,waveA);
%{ 
use following code if you want to calculate iso12233
Ascene = sceneCreate('uniform equal energy',nPixels,wave);
scene = sceneCreate('slantededge',62,2.6,[],wave);
Ascene = sceneAdjustIlluminant(Ascene,sceneIlluminant);
scene = sceneAdjustIlluminant(scene,sceneIlluminant);
%scene = sceneAdd(scene,Ascene,'add'); %put the two together
%scene = Ascene; %just go without slantededge
%}
% Calculate the fluorescence for this illuminant
illuminant = sceneGet(scene,'illuminant photons');
%Choose your Fluorophore here:
fName  = fullfile(fiToolboxRootPath,'data','Monici','Porphyrins.mat');
fl  = fiReadFluorophore(fName,'wave',waveA);
sz = sceneGet(scene,'size');
fLuminance = 2;
%% Instantiate a Camera
sensor = sensorCreate('bayer (gbrg)');  % just the common Bayer one
% Go with Joyce's Defaults
voltageSwing   = 1.15;  % Volts
wellCapacity   = 9000;  % Electrons
conversiongain = voltageSwing/wellCapacity;   
fillfactor     = 0.9;       % A fraction of the pixel area
pixelSize      = 2.2*1e-6;   % Meters
darkvoltage    = 1e-005;     % Volts/sec
readnoise      = 0.00096;    % Volts
sensorSet(sensor,'pixel size',[pixelSize pixelSize]*fillfactor);
sensor = sensorSet(sensor,'pixel size constant fill factor',[pixelSize pixelSize]);
sensor = sensorSet(sensor,'pixel conversion gain',conversiongain);
sensor = sensorSet(sensor,'pixel voltage swing',voltageSwing);
sensor = sensorSet(sensor,'pixel dark voltage',darkvoltage);
sensor = sensorSet(sensor,'pixel read noise volts',readnoise);
%Now make the lens
oi = oiCreate;
oi = oiSet(oi,'optics fnumber',4);
oi = oiSet(oi,'optics offaxis','cos4th');
oi = oiSet(oi,'optics focal length',3e-3);
exposureDuration = 0.5; % comment this part out if you're gonna auto-expose
dsnu =  0.0010;           % Volts (dark signal non-uniformity)
prnu = 0.2218;            % Percent (ranging between 0 and 100) photodetector response non-uniformity
analogGain   = 1;         % Used to adjust ISO speed
analogOffset = 0;         % Used to account for sensor black level
rows = 466;               % number of pixels in a row
cols = 642;               % number of pixels in a column
sensor = sensorSet(sensor,'exposuretime',exposureDuration); % comment in if you want to set exposure time
sensorSet(sensor,'autoExposure',0);  % set to zero if you ain't auto-exposing
sensor = sensorSet(sensor,'rows',rows);
sensor = sensorSet(sensor,'cols',cols);
sensor = sensorSet(sensor,'dsnu level',dsnu);  
sensor = sensorSet(sensor,'prnu level',prnu); 
sensor = sensorSet(sensor,'analog Gain',analogGain);     
sensor = sensorSet(sensor,'analog Offset',analogOffset); 
%% Now watch the magic
tic;
for i = 13 %1:101 choose indices of the lowpass cutoff that you want
for j = 12 %1:31 choose indices of the highpass cutoff that you want
l_co = lowpass_cutoff(i);
h_co = highpass_cutoff(j);
%lowpass filter, modeled by a gaussian
cutoff = l_co;
sigma = 20;
x = linspace(380,780,101);
[minValue,closestIndex] = min(abs(cutoff-x));
closestValue = x(closestIndex);
filter = normpdf(x,cutoff,sigma) + 0.001;
filter(1:closestIndex) = max(filter);
sceneIlluminantNew = sceneIlluminant.*filter;
%plot(linspace(380,780,101),sceneIlluminant);
scene = sceneAdjustIlluminant(scene,sceneIlluminantNew);
%sceneWindow(scene);
slope = 2.6;
pattern = imageSlantedEdge(sz-1,slope);
fScene = fiSceneCreate(fl,pattern,sceneIlluminantNew);
fScene = sceneAdjustLuminance(fScene,fLuminance);
%sceneWindow(fScene);
%% Combine the scenes
% Add the radiance due to reflectance and the radiance due to fluorescence
cScene = fiSceneAddFluorescence(scene, fScene );
cScene = sceneSet(cScene,'name','Combined');
%sceneWindow(cScene);
oi = oiCompute(cScene,oi);
%% color filters and highpass now
% Change this to be the Sony sensor
wave = sensorGet(sensor,'wave');
fullFileName = fullfile(isetRootPath,'data','sensor','colorfilters','nikon','NikonD100.mat');
[data3,filterNames] = ieReadColorFilter(wave,fullFileName); 
%generate a highpass filter using a gaussian
cutoff = h_co;
sigma = 10;
freq = linspace(wave(1),wave(end),length(wave));
[minValue,closestIndex] = min(abs(cutoff-freq));
closestValue = freq(closestIndex);
longpass_filter = normpdf(freq,cutoff,sigma);
longpass_filter(closestIndex:length(longpass_filter)) = max(longpass_filter);
data3(:,1) = data3(:,1).*longpass_filter';
data3(:,2) = data3(:,2).*longpass_filter';
data3(:,3) = data3(:,3).*longpass_filter';
sensor = sensorSet(sensor,'filter spectra',data3);
sensor = sensorSet(sensor,'filter names',filterNames);
sensor = sensorSet(sensor,'Name','Camera-Simulation');
%% compute sensor image
sensor = sensorCompute(sensor,oi);
%sensorWindow(sensor);
%% Image process
ip = ipCreate;
% To see it clearly, you can do this
ip = ipCompute(ip,sensor);
ipWindow(ip); %uncomment if you wanna see the picture
%get a metric for contrast
dispimg = ip.data.result;
dispimg = dispimg(113:353,200:442,:); %crop the data
d_metric = mean2(dispimg(53:70,58:64,1)) + mean2(dispimg(53:70,58:64,2)) + mean2(dispimg(53:70,58:64,3));
d_metric = d_metric - mean2(dispimg(201:215,224:229,1)) - mean2(dispimg(201:215,224:229,2)) - mean2(dispimg(201:215,224:229,3));
contrast_metric(i,j) = d_metric; %uncomment to measure
%buffer(b) = d_metric %uncomment when sweeping
end
end
toc;
end
%% generate the heat map

%uncomment if you wanna see the heat map
%pcolor(highpass_cutoff,lowpass_cutoff,contrast_metric);
%colorbar;
%xlabel('Highpass Cutoff Frequency');
%ylabel('Lowpass Cutoff Frequency');

%uncomment if you are sweeping across light gaussians
%plot(light_mus,buffer);

%save your parameters if you wanna use them alter
%saveas(gcf,'contrast2collagen.png');
%save('cmetric_collagen.mat','contrast_metric');