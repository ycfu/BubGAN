clc
clearvars -except BubImage
close all
warning off
%Set the plot parameters**********************************
nF=1;
set(0,'DefaultFigureUnits','pixels','DefaultFigurePosition',[0 0 round(nF*800) round(nF*800)])
set(0,'DefaultFigureColor',[1 1 1])
set(0,'DefaultAxesUnits','normalized','DefaultAxesPosition',[0.13 0.13 0.75 0.8])
%set(0,'DefaultAxesXTickMode','manual','DefaultAxesYTickMode','manual')
set(0,'DefaultAxesTickLength',[0.02 0.02])
set(0,'DefaultAxesXMinorTick','off','DefaultAxesYMinorTick','off')
set(0,'DefaultAxesLineWidth',ceil(1.5),'DefaultAxesFontName','Arial',...
    'DefaultAxesFontSize',ceil(30),'DefaultAxesBox','on')
set(0,'DefaultLineLineWidth',ceil(1.5),'DefaultLineMarkerSize',ceil(8))
set(0,'DefaulttextFontName','Arial','DefaulttextFontSize',26)
%Set the plot parameters**********************************

rng(10)  % for reproducibility
addpath ./BubFunctions
addpath ./ColorMap

InputFolder=['.\MillionBubbleDataBase\'];

if ~exist('BubImage','var')
    load([InputFolder,'MillionBub.mat'])  % load million bubble image database
end
SelectBubRange=1000000;                     % choose from 1:1000000 images
load([InputFolder,'Synthetic_BubInfo.mat']) % Corresponding bubble information

ImageGenNum=1; % Generating image number
SaveFlag=1;    % Do not save data

%% Parameter setup for bubbly flow generation
%% Setup interested simulation domain
ParaImg.Depth=10; % mm
ParaImg.Width=30; % mm
ParaImg.Height=60; % mm
ParaImg.VoidFraction=0.1;  %[-] The void fraction is calculated by estimating third axis along camera optical axis is same to the short axis in 2D projection
ParaImg.pixtomm=0.05;   % mm/pix  Image Resolution
ParaImg.bckgrdint=183;
ParaImg.TotVol=ParaImg.Depth*ParaImg.Width*ParaImg.Height;   % mm^3

%% Bubble size distribution
ParaImg.BubRefSize=3.5;  % Bubble mean diameter
ParaImg.BubSizeMode='uniform'; % uniform or Gaussian
ParaImg.BubDev=1;  % If uniform distribution,  bubble diameter range from [BubRefSize-BubDev, BubRefSize+BubDev]
% If Gaussian distribution, BubDev serve as standard deviation

%% Bubble Location distribution along x direction
ParaImg.DistX.Flag=1;
ParaImg.DistX.Loc=[1:1:30];
ParaImg.DistX.Counts(1:30)=10;  % uniform density. The value is relative value.
%
% ParaImg.DistX.Counts(1:5)=10;  %This is non-uniform bubble density 
% ParaImg.DistX.Counts(6:25)=4;
% ParaImg.DistX.Counts(26:30)=10;

%% Generate Image Canvas
% The Image Canvas can also be loaded from existed background
ImageCanvas=double(ones(round(ParaImg.Height/ParaImg.pixtomm),round(ParaImg.Width/ParaImg.pixtomm)))*ParaImg.bckgrdint;
ImageCanvas=ImageCanvas+2*randn(size(ImageCanvas)); % Add noise
ParaImg.ImageCanvas=uint8(ImageCanvas);

%% Image Generation
for k=1:ImageGenNum
    disp(['Generating the image of ', num2str(k,'%3.f')])
    [BubList]=BubDistribution( ParaImg,[ParaImg.Height, ParaImg.Width, ParaImg.Depth]); % generated bubble parameters
    [ImageCanvas_Painted,ImgLabel]=DrawCanvas(BubList,ParaImg,BubImage,BubInfo,SelectBubRange); %draw bubbles on canvas
    if SaveFlag
    %% Save Generated images and bubble label information
    imwrite(ImageCanvas_Painted, sprintf('Image_%03.f.tif', k));
    save( sprintf('Label_%03.f.mat',k),'ImgLabel')
    end
end

%% Display Painted Image
figure,
imshow(ImageCanvas_Painted)

%% Area Label information plot

figure,
patch([0 0 size(ImageCanvas_Painted,2) size(ImageCanvas_Painted,2) 0],[0 size(ImageCanvas_Painted,1) size(ImageCanvas_Painted,1) 0 0],[0 0 0],'FaceAlpha',0.2)
hold on
for k=1:length(ImgLabel)
    boundary=ImgLabel(k).boundary;
    index1=find(boundary(:,2)<0);
    index2=find(boundary(:,2)>size(ImageCanvas_Painted,2));
    index3=find(boundary(:,1)<0);
    index4=find(boundary(:,1)>size(ImageCanvas_Painted,1));
    boundary([index1 index2 index3 index4],:)=[];
    
    area = polyarea(boundary(:,1),boundary(:,2));
    area = area*ImgLabel(k).resolution^2;
    patch(boundary(:,2),boundary(:,1),area,'FaceAlpha',1)
    plot(ImgLabel(k).xx,ImgLabel(k).yy,'r.')
end
plot([0 0 size(ImageCanvas_Painted,2) size(ImageCanvas_Painted,2) 0],[0 size(ImageCanvas_Painted,1) size(ImageCanvas_Painted,1) 0 0],'k','LineWidth',1)
axis equal
axis off
set(gca, 'YDir','reverse')
colormap(brewermap([],'*Spectral'))
c=colorbar
ylabel(c,'Area $[\mathrm{mm}^2]$','Interpreter','Latex')
if SaveFlag, saveas(gcf,'Bubbly_Labeled_Area.tif'),end
%% Rotation angle Label information plot
figure,
patch([0 0 size(ImageCanvas_Painted,2) size(ImageCanvas_Painted,2) 0],[0 size(ImageCanvas_Painted,1) size(ImageCanvas_Painted,1) 0 0],[0 0 0],'FaceAlpha',0.2)
hold on
for k=1:length(ImgLabel)
    boundary=ImgLabel(k).boundary;
    index1=find(boundary(:,2)<0);
    index2=find(boundary(:,2)>size(ImageCanvas_Painted,2));
    index3=find(boundary(:,1)<0);
    index4=find(boundary(:,1)>size(ImageCanvas_Painted,1));
    boundary([index1 index2 index3 index4],:)=[];
    patch(boundary(:,2),boundary(:,1),ImgLabel(k).BubInfo.phi,'FaceAlpha',1)
    plot(ImgLabel(k).xx,ImgLabel(k).yy,'r.')
end
plot([0 0 size(ImageCanvas_Painted,2) size(ImageCanvas_Painted,2) 0],[0 size(ImageCanvas_Painted,1) size(ImageCanvas_Painted,1) 0 0],'k','LineWidth',1)
axis equal
axis off
set(gca, 'YDir','reverse')
colormap(brewermap([],'PiYG'))
c=colorbar;
caxis([-pi/2 pi/2])
ylabel(c,'Rotation angle $\varphi$ [rad] ','Interpreter','Latex')
if SaveFlag, saveas(gcf,'Bubbly_Labeled_phi.tif'),end

%% Aspect Ratio information plot
figure,
patch([0 0 size(ImageCanvas_Painted,2) size(ImageCanvas_Painted,2) 0],[0 size(ImageCanvas_Painted,1) size(ImageCanvas_Painted,1) 0 0],[0 0 0],'FaceAlpha',0.2)
hold on
for k=1:length(ImgLabel)
    boundary=ImgLabel(k).boundary;
    index1=find(boundary(:,2)<0);
    index2=find(boundary(:,2)>size(ImageCanvas_Painted,2));
    index3=find(boundary(:,1)<0);
    index4=find(boundary(:,1)>size(ImageCanvas_Painted,1));
    boundary([index1 index2 index3 index4],:)=[];
    patch(boundary(:,2),boundary(:,1),ImgLabel(k).BubInfo.E,'FaceAlpha',1)
    plot(ImgLabel(k).xx,ImgLabel(k).yy,'r.')
end
plot([0 0 size(ImageCanvas_Painted,2) size(ImageCanvas_Painted,2) 0],[0 size(ImageCanvas_Painted,1) size(ImageCanvas_Painted,1) 0 0],'k','LineWidth',1)
axis equal
axis off
set(gca, 'YDir','reverse')
colormap(brewermap([],'*YlOrBr'))
c=colorbar;
caxis([0.5 1])
ylabel(c,'Aspect ratio $E$ [-] ','Interpreter','Latex')
if SaveFlag, saveas(gcf,'Bubbly_Labeled_E.tif'),end

%% Bubble Intensity information plot
XCenter=[ImgLabel.xx]';
YCenter=[ImgLabel.yy]';
Center=[XCenter YCenter];
C=zeros(size(ImageCanvas_Painted));
Sigma = [100 0; 0 100];
x1 =1:1:600; x2 =1:1:1200;
[X1,X2] = meshgrid(x1,x2);
for k=1:size(XCenter)
mu = [XCenter(k) YCenter(k)];
F = mvnpdf([X1(:) X2(:)],mu,Sigma);
F = reshape(F,length(x2),length(x1));
C=C+F;
end
figure,
imagesc(x1,x2,C);
axis equal
axis off
if SaveFlag, saveas(gcf,'Bubbly_Labeled_Intensity.tif'),end


