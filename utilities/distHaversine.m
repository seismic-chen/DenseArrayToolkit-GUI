function [azeq,azst,delta,dist] = distHaversine(eqlat,eqlon,stlat,stlon)
% modified from Dr. Jeff Gu's finddelta code
% calculate the great arc distance from the event to the station
% function [azeq,azst,delta,dist] = distHaversine(eqlat,eqlon,stlat,stlon)
rad = 0.0174533;
reprad = 57.295779;
geoco = 0.993277;
sq302 = 0.8660254;
if abs(stlat) >= 90
stlat = sign(stlat)*89.98;
end
stl = atan(geoco*tan(stlat*rad));
slon = stlon*rad;
as = sin(stl);
bs = cos(stl);
cs = cos(slon);
ds = sin(slon);

elat = atan(geoco*tan(eqlat*rad));
elon = eqlon*rad;
a = sin(elat);
b = cos(elat);
c = cos(elon);
d = sin(elon);
f1 = 1.5*a^2-0.5;
f2 = 2.0*a*b*sq302;
f3 = b^2*sq302;

cdif = c*cs+d*ds;
sdif = c*ds-d*cs;
cdel = a*as+b*bs*cdif;
if abs(cdel) > 1
cdel = 1.0*sign(cdel);
end
azeq = reprad*atan2(bs*sdif,(as-a*cdel)/b);
azst = reprad*atan2(-b*sdif,(a-as*cdel)/bs);
% degree
delta = reprad*acos(cdel);
% degree to distance
dist=delta*6371.0*3.141592653579/180.0;
if azeq < 0
azeq = 360.0+azeq;
end
if azst < 0
azst = 360.0+azst;
end
%sprintf('station azimuth: %10.6f',azst)
%sprintf('station azimuth: %10.6f',azeq)