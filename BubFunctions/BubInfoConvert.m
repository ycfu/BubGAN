function [BubInfoNorm] = BubInfoConvert()
%BUBINFOCONVERT Summary of this function goes here
%   Output the normailized BubInfo

BubInfo=load('Run1_C1_Info_Org_BubInfo_V2.mat');
BubInfo=BubInfo.BubInfo;
BubInfo=table2array(BubInfo);
% 1.rotation angle phi (-pi/2, pi/2]
% 2. Aspect ratio E (0, 1], 3. Edge ratio (0,1]
% 4.Circularity Psi
BubInfo=BubInfo(:,[3 6 8 9])';
for k=1:size(BubInfo,1)
%     BubInfo(k,:)=rescale(BubInfo(k,:),-1,1);
 BubInfoNorm(k,:)=( BubInfo(k,:)-mean( BubInfo(k,:)))/std( BubInfo(k,:));
end


end

