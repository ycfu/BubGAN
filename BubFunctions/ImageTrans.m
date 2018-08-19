function [I1]=ImageTrans(I,option)

if strcmp(option,'peel')
    if mod(size(I,1),2)==0
        R=size(I,1); C=size(I,2);
        linearInd=[];
        linearIndPart=[];
        for layer=1:size(I,1)/2
            R1=layer; R2=size(I,2)-layer+1;
            C1=layer; C2=size(I,1)-layer+1;
            Rvec=[R1:1:R2];
            Cvec=[C1:1:C2];
            linearIndA = sub2ind([R C], Rvec(1)*ones(size(Cvec)), Cvec);
            linearIndB = sub2ind([R C], Rvec(2:end-1),  Cvec(end)*ones(size(Rvec(2:end-1))));
            linearIndC = sub2ind([R C], Rvec(end)*ones(size(Cvec)), flip(Cvec) );
            linearIndD = sub2ind([R C], flip(Rvec(2:end-1)), Cvec(1)*ones(size(Rvec(2:end-1))));
            linearInd=[linearInd linearIndA linearIndB linearIndC linearIndD];
        end
        I1=I(linearInd);
        I1=reshape(I1,size(I));
        for k=2:2:size(I,1)
            I1(:,k)=flip(I1(:,k));
        end
    end
end


if strcmp(option,'wrap')
    if mod(size(I,1),2)==0
        for k=2:2:size(I,1)
            I(:,k)=flip(I(:,k));
        end
        R=size(I,1); C=size(I,2);
        linearInd=[];
        linearIndPart=[];
        for layer=1:size(I,1)/2
            R1=layer; R2=size(I,2)-layer+1;
            C1=layer; C2=size(I,1)-layer+1;
            Rvec=[R1:1:R2];
            Cvec=[C1:1:C2];
            linearIndA = sub2ind([R C], Rvec(1)*ones(size(Cvec)), Cvec);
            linearIndB = sub2ind([R C], Rvec(2:end-1),  Cvec(end)*ones(size(Rvec(2:end-1))));
            linearIndC = sub2ind([R C], Rvec(end)*ones(size(Cvec)), flip(Cvec) );
            linearIndD = sub2ind([R C], flip(Rvec(2:end-1)), Cvec(1)*ones(size(Rvec(2:end-1))));
            linearInd=[linearInd linearIndA linearIndB linearIndC linearIndD];
        end
        I1(linearInd)=I;
        I1=reshape(I1,size(I));
    end
end

end
