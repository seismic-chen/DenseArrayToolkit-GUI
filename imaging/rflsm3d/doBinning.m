function [ dout,mask] = doBinning(din,rx,ry, x, y,nx,ny, dx, dy)

% input:
%       din: nt*ntraces matrix
%       rx: x stations coordinates
%       ry: y stations coordinates
%       x: x regular grid coordinates
%       y: y regular grid coordinates
%       nx: points of x axis
%       ny: points of y axis
%       dx: spacing of x
%       dy: spacing of y
%
% note: Stations' coordinates MUST align with each column of din

[n1,n2]=size(din);

dout=zeros(n1,nx,ny);
mask=ones(nx,ny);

for iy=1:ny
    for ix=1:nx
        index=find(rx>=x(ix)-dx/2 & rx<=x(ix)+dx/2 & ry>=y(iy)-dy/2 & ry<=y(iy)+dy/2);
        n=length(index);
        if n==1
            dout(:,ix,iy)=din(:,index);
        elseif n==0
            mask(ix,iy)=0;
            dout(:,ix,iy)=zeros(n1,1);
        end

        if n>=2
            if rx(index(1))==x(ix) && ry(index(1))==y(iy)
                dout(:,ix,iy)=din(:,index(1));
            else
                % t1=sqrt((x(index(1))-xout(ix))^2+(y(index(1))-yout(iy))^2);
                % t2=sqrt((x(index(2))-xout(ix))^2+(y(index(2))-yout(iy))^2);
                % dout(:,ix,iy)=(t1*din(:,index(2))+t2*din(:,index(1)))/(t1+t2);
                [dist,idx]=sort(sqrt((rx(index)-x(ix)).^2+(ry(index)-y(iy)).^2));
                t1 = dist(1);
                t2 = dist(2);

                dout(:,ix,iy)=(t1*din(:,index(idx(2)))+t2*din(:,index(idx(1))))/(t1+t2);
                %                     dout(:,ix)=(t1*din(:,index(2))+t2*din(:,index(1)))/(t1+t2);
            end
        end
    end
end

if nargout==4
    mask=ones(n1,1)*reshape(mask,1,nx*ny);
    mask=reshape(mask,n1,nx,ny);
end

return