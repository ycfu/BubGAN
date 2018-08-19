function  [BubList]=BubDistribution( ParaImg, Dim)
Height=Dim(1);
Width=Dim(2);
Depth=Dim(3);
ChannelCapacity=Height*Width*Depth;


TotalVolume=0;
BubRefSize=ParaImg.BubRefSize;
BubDev=ParaImg.BubDev;
i=0;

while TotalVolume<ChannelCapacity*ParaImg.VoidFraction
    i=i+1;
    if strcmp(ParaImg.BubSizeMode,'uniform')
        a=BubRefSize/2+BubDev/2*(-1+2*rand);
        b=BubRefSize/2+BubDev/2*(-1+2*rand);
    elseif strcmp(ParaImg.BubSizeMode,'Gaussian')
        a=normrnd(BubRefSize/2,BubDev/2);
        b=normrnd(BubRefSize/2,BubDev/2);
        while a<0 | b<0   % reject negative semi-axis and redo
            a=normrnd(BubRefSize/2,BubDev/2);
            b=normrnd(BubRefSize/2,BubDev/2);
        end
    end
    c=min(a,b);
    maxbound=max([a b c]);
    if b>a
       temp=a;
       a=b;
       b=temp;
    end
    phi=pi/2*(0.5-rand);
    TotalVolume=TotalVolume+4/3*pi*a*b*c;
    BubList(i).a=a;
    BubList(i).b=b;
    BubList(i).c=c;
    BubList(i).phi=phi;
end

temp= struct2table(BubList);
temp=table2array(temp);
a_avg=mean(temp(:,1));
b_avg=mean(temp(:,2));
semi_avg=(a_avg+b_avg)/2;
semi_axes=temp(:,1:3);


if ParaImg.DistX.Flag==0
    xlist= Width*rand(1,length(BubList));
else
    xlist=SamplefromCounts(ParaImg,length(BubList));
end
ylist= Height*rand(1,length(BubList));
zlist= Depth*rand(1,length(BubList));
loc=[xlist' ylist' zlist'];

for k=1:length(BubList)   % Check if bubble physically touched or overlap together.
    if k>1
        D1=0;
        D2=0;
        while min(D1-D2)<0
            D1=pdist2(loc(1:k-1,:),semi_axes(k,:),'euclidean');
            D2=max(semi_axes(k-1,:),max(semi_axes(k,:)))*2;
            if min(D1-D2)<0
                loc(k,2)=semi_avg+(Height-2*semi_avg)*rand;
                loc(k,3)=semi_avg+(Depth-2*semi_avg)*rand;
            end
        end
    end
end

for k=1:length(BubList)
    BubList(k).xx=loc(k,1);
    BubList(k).yy=loc(k,2);
    BubList(k).zz=loc(k,3);
end
end

%% Sample Location of bubble number density distribtion
function [Loc]=SamplefromCounts(ParaImg,nSample)

LocinPix=[ParaImg.pixtomm/2:ParaImg.pixtomm:ParaImg.Width-ParaImg.pixtomm/2];
CountsinPix=interp1(ParaImg.DistX.Loc,ParaImg.DistX.Counts,LocinPix,'linear','extrap');
CountsinPix=CountsinPix/max(1e-20,min(CountsinPix));
CountsinPixcum=cumsum(CountsinPix);
CP=floor(rand(1,nSample)*CountsinPixcum(end));
for k=1:nSample
    Loc(k)=LocinPix(find(CP(k)<=CountsinPixcum,1));
end
end
