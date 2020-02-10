function [nugap,freq,distnu,w]=richgap(L1, L2, p1, p2)

L1=minreal(L1);
L2=minreal(L2);


%Important: here it defines the frequencies of interest to compute the
%nu-gap
w=(p2:0.001:p1);

%Now it decompose the transfer function in immmaginary and real part in the
%nyquist plane
[re,im]=nyquist(L1,w);
re=squeeze(re);im=squeeze(im);
[re1,im1]=nyquist(L2,w);
re1=squeeze(re1);im1=squeeze(im1);

%Here it defines the riemann sphere on which it will project the nyquist
%plots
[X,Y,Z]=sphere(50);
X=X/2;Y=Y/2;Z=0.5+Z/2;


%It computes the projections on the sphere of the transfer functions
[x,y,z]=RiemannProject(re,im);
[x1,y1,z1]=RiemannProject(re1,im1);


%Finally it computes the maximum between the distances

[nugap,I,distnu]=maxdistance([x,y,z],[x1,y1,z1]);
freq=w(I);



%take real and imaginary parts of a nyquist contour and project them onto a
%reimann sphere
function [x,y,z]=RiemannProject(re,im)

%now projection is much easier in polars:
[theta,r]=cart2pol(re,im);
%define the z point of the projection
z=(r.^2)./(1+r.^2);
%the relevant polar r is then given by
rproj=r.*(1-z);
%now shift into cartesian again
[x,y]=pol2cart(theta,rproj);

end


function [gap,II,dist]=maxdistance(a,b)

for i=1:size(a,1)    
    dist(i)=sqrt((a(i,1)-b(i,1))^2+(a(i,2)-b(i,2))^2+(a(i,3)-b(i,3))^2);
end;
    
[gap, II]=max(dist);


end


end