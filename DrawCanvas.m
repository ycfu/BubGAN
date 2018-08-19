function [ImageCanvas,ImgLabel]=DrawCanvas(BubList,ParaImg,BubImage,BubInfo,SelectBubRange)
ImageCanvas=ParaImg.ImageCanvas;

for k=1:size(BubList,2)
    disp(['      Generating bubble ', num2str(k),' of ', num2str(size(BubList,2))])
    RandBubNum=randi([1 SelectBubRange]);
    patchimage = BubImage(:,:,RandBubNum);
    [patchBubInfo patchboundary]=BubInfoExt(patchimage);
    [patchimage patchBubInfo patchboundary mask mask1 mask2 ratio]=patchscale(patchimage, patchBubInfo, patchboundary, BubList(k).a/ParaImg.pixtomm);
    
    % center (xxg yyg) and boundary (xg yg) in global canvas with pixel unit
    xg=patchboundary(:,2)+BubList(k).xx/ParaImg.pixtomm-patchBubInfo.xx;
    yg=patchboundary(:,1)+BubList(k).yy/ParaImg.pixtomm-patchBubInfo.yy;
    xxg=BubList(k).xx/ParaImg.pixtomm;
    yyg=BubList(k).yy/ParaImg.pixtomm;
    
    % apply boundary physical restriction on x direction
    [xxg, yyg, xg,yg]=phyrestriction(ParaImg.ImageCanvas,xxg, yyg, xg,yg, 'x');
    %[xxg, yyg, xg,yg]=phyrestriction(ParaImg.ImageCanvas,xxg, yyg, xg,yg, 'y');
    
    locx=round(xxg-patchBubInfo.xx)+1;
    locx1=locx+size(patchimage,1)-1;
    locy=round(yyg-patchBubInfo.yy)+1;
    locy1=locy+size(patchimage,1)-1;
    patchloc=[locx locx1 locy locy1];
    
    
    [ImageCanvas, patchimage, mask, mask1, mask2, patchloc]=patchcrop(ImageCanvas, patchimage, mask, mask1, mask2, patchloc);
    ImgLabel(k).BubInfo=patchBubInfo;
    ImgLabel(k).boundary=[yg xg];
    ImgLabel(k).xx=xxg;
    ImgLabel(k).yy=yyg;
    ImgLabel(k).resolution=ParaImg.pixtomm;
    ImgLabel(k).patchimage=patchimage;
end
end






%% Fine physical restriction
function [xcenter, ycenter, x,y]=phyrestriction(ImageCanvas,xcenter, ycenter, x,y, direction)
if strcmp(direction,'x')==1
    if max(x(:))>size(ImageCanvas,2)
        xcenter=xcenter- (max(x(:))-size(ImageCanvas,2));
        x=x-(max(x(:))-size(ImageCanvas,2));
    end
    if min(x(:))<0
        xcenter=xcenter+1-min(x(:));
        x=x+1-min(x(:));
    end
end
if strcmp(direction,'y')==1
    if max(y(:))>size(ImageCanvas,1)
        ycenter=ycenter- (max(y(:))-size(ImageCanvas,1));
        y=y-(max(y(:))-size(ImageCanvas,1));
    end
    if min(y(:))<0
        ycenter=ycenter+1- min(y(:));
        y=y+1-min(y(:));
    end
end
end


%% patchscale scale the patch to given size
function [patchimage patchBubInfo patchboundary mask mask1 mask2 ratio]=patchscale(patchimage, patchBubInfo, patchboundary, Bubsize)
% Bubsize is the given size for bubble generation.  The 64 by 64 pixel
% bubble image from million-bubble database should be scaled to the
% given size.
% Bubsize should be given in pixel.

ratio=Bubsize/patchBubInfo.aa;
patchpix=size(patchimage,1)*ratio;
patchimage=imresize(patchimage, [patchpix patchpix]);
patchBubInfo.aa = patchBubInfo.aa*ratio;
patchBubInfo.bb = patchBubInfo.bb*ratio;
patchBubInfo.xx = patchBubInfo.xx*ratio;
patchBubInfo.yy = patchBubInfo.yy*ratio;
patchboundary=double(patchboundary*ratio);
mask=poly2mask(patchboundary(:,2),patchboundary(:,1),size(patchimage,1),size(patchimage,2));
mask1=bwmorph(mask,'thicken',1);
mask2=bwmorph(mask,'thicken',2)- bwmorph(mask,'thicken',0);
end

%% patchcrop: crops the patch parts which is outside the ImageCanvas and then paint it
function [ImageCanvas, patchimage, mask, mask1, mask2, patchloc]=patchcrop(ImageCanvas, patchimage, mask, mask1, mask2, patchloc)
dx=0; dx1=0; dy=0; dy1=0;
locx=patchloc(1);
locx1=patchloc(2);
locy=patchloc(3);
locy1=patchloc(4);
if locx<=0
    dx=1-locx;
    locx=locx+dx;
    mask1=mask1(:,1+dx:end);
    mask2=mask2(:,1+dx:end);
    patchimage=patchimage(:,1+dx:end);
end
if locx1>size(ImageCanvas,2)
    dx1=locx1-size(ImageCanvas,2);
    locx1=locx1-dx1;
    mask1=mask1(:,1:end-dx1);
    mask2=mask2(:,1:end-dx1);
    patchimage=patchimage(:,1:end-dx1);
end
if locy<=0
    dy=1-locy;
    locy=locy+dy;
    mask1=mask1(1+dy:end,:);
    mask2=mask2(1+dy:end,:);
    patchimage=patchimage(1+dy:end,:);
end
if locy1>size(ImageCanvas,1)
    dy1=locy1-size(ImageCanvas,1);
    locy1=locy1-dy1;
    mask1=mask1(1:end-dy1,:);
    mask2=mask2(1:end-dy1,:);
    patchimage=patchimage(1:end-dy1,:);
end

maskcanvas1=boolean(zeros(size(ImageCanvas)));
maskcanvas1(locy:locy1,locx:locx1)=mask1;
maskcanvas2=boolean(zeros(size(ImageCanvas)));
maskcanvas2(locy:locy1,locx:locx1)=mask2;
patchimage=patchimage+uint8(double(mode(ImageCanvas(maskcanvas1)))-double(mode(patchimage(:)))); % adjust bubble intensity according to surrounding background(ImageCanvas)
ImageCanvas(maskcanvas1)=patchimage(mask1);
ImageCanvas=regionfill(ImageCanvas,maskcanvas2);

end



