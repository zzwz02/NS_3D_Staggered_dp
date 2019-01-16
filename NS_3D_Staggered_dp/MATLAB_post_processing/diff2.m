function B = diff2(A,k,bc)

if size(A,1)==1, A = A'; end
[nx,ny,nz]=size(A);
if (bc==1)
    if (k==1)
        A1=zeros(nx+1,ny,nz);
        A1(2:end,:,:) = A;
        A1(1,:,:) = A(end-1,:,:);
    elseif (k==2)
        A1=zeros(nx,ny+1,nz);
        A1(:,2:end,:) = A;
        A1(:,1,:) = A(:,end-1,:);
    elseif (k==3)
        A1=zeros(nx,ny,nz+1);
        A1(:,:,2:end) = A;
        A1(:,:,1) = A(:,:,end-1);
    end
    A=A1;
end

if (k==1)
    B = A(1:end-2,:,:)-2*A(2:end-1,:,:)+A(3:end,:,:);
elseif (k==2)
    B = A(:,1:end-2,:)-2*A(:,2:end-1,:)+A(:,3:end,:);
elseif (k==3)
    B = A(:,:,1:end-2)-2*A(:,:,2:end-1)+A(:,:,3:end);
end
if size(A,2)==1, B = B'; end