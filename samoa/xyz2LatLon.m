function [Lon,Lat] = xyz2LatLon(X,Y)

%reshape
x = reshape(X,[],1);
y = reshape(Y,[],1);

%utmZone
zone='27 W';
A = repmat(zone,length(x),1);

%Conversion
[Lat, Lon]=utm2deg(x,y,A); 

%reshape pt2
Lat = reshape(Lat, sqrt(length(Lat)), []);
Lon = reshape(Lon, sqrt(length(Lon)), []);

end