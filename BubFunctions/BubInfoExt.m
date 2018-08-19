function [BubInfo boundary]= BubInfoExt(I)


%% Extract bubble boundary from the image patch
if size(I,3)==3
    I=rgb2gray(I);
end

bw = im2bw(I);
bw=1-bw;
bwb= bwmorph(bw,'bridge',Inf);
bwb= bwmorph(bwb,'bridge',Inf);
bwb_withhole=bwb;
edgepix=length(find(bwb_withhole==1));
bwb=imfill(bwb,'holes');
filledpix=length(find(bwb==1));
bwb=bwmorph(bwb,'thicken',5);
bwb=activecontour(I, bwb, 10, 'edge');

fgind=(bwb_withhole==1);
bgind=(bwb_withhole==0);

FgInt=mean(I(fgind));
BgInt=mean(I(bgind));




[B]=bwboundaries(bwb,8,'noholes');      % choose the longest boundary if there are two or more
if length(B)>1
  len=cellfun(@(x) size(x,1),B);
  boundary=B{find(len==max(len))}; 
else
  boundary=B{1};
end

tempx=boundary(:,2);
tempy=boundary(:,1);
dx=abs(tempx-circshift(tempx,1));
dy=abs(tempy-circshift(tempy,1));

boundlen=sum(sqrt(dx.^2+dy.^2));  
Psi=4*pi*filledpix/boundlen^2;

%% Fit the boundary with ellipse to extract orientatin and semi-axes

[xx yy aa bb para ellipse_t] = ellipsefitting( boundary(:,2),boundary(:,1) );

phi=ellipse_t.phi;    % convert the rotation angle to [-pi/2 pi/2], which it the angle between long semi-axis and horizontal direction
if aa<bb
    temp=aa;
    aa=bb;
    bb=temp;
    phi=phi+pi/2;
    if phi>pi/2
        phi=phi-pi;
    end
end

[Gmag, Gdir] = imgradient(I);
Gmag=Gmag/8; % normalize the intenstiy gradient to level/pix.

idx=sub2ind(size(Gmag),boundary(:,2),boundary(:,1));
indices=sub2ind(size(Gmag),boundary(:,1),boundary(:,2));
for k=1:length(boundary(:,1))
    blockszie=4;
    boundarygradient(k)=max(max(  Gmag( max(1,boundary(k,1)-blockszie):min(size(Gmag,1),boundary(k,1)+blockszie), ...
        max(1,boundary(k,2)-blockszie):min(size(Gmag,2),boundary(k,2)+blockszie)) ));
end
EdgeGrad=mean(boundarygradient);

%%  Assign to BubInfo var.

BubInfo=table(single(aa),single(bb), single(phi), single(xx), single(yy), single(bb/aa),...
    single(EdgeGrad), single(edgepix/filledpix),single(Psi),...
    single(FgInt),single(BgInt),...
    'VariableNames',{'aa','bb','phi','xx','yy','E','EdgeGrad','ER','Psi','FgInt','BgInt'});


end



