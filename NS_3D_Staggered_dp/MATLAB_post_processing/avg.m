function B = avg(A,k,bc)

if size(A,1)==1, A = A'; end
if (k==1)
    B = (A(2:end,:,:)+A(1:end-1,:,:))/2;
elseif (k==2)
    B = (A(:,2:end,:)+A(:,1:end-1,:))/2;
elseif (k==3)
    B = (A(:,:,2:end)+A(:,:,1:end-1))/2;
end

[nx,ny,nz]=size(B);

if (bc==1)
    if (k==1)
        B1=zeros(nx+1,ny,nz);
        B1(2:end,:,:) = B;
        B1(1,:,:) = B(end,:,:);
    elseif (k==2)
        B1=zeros(nx,ny+1,nz);
        B1(:,2:end,:) = B;
        B1(:,1,:) = B(:,end,:);
    elseif (k==3)
        B1=zeros(nx,ny,nz+1);
        B1(:,:,2:end) = B;
        B1(:,:,1) = B(:,:,end);
    end
    B=B1;
end
if size(A,2)==1, B = B'; end
