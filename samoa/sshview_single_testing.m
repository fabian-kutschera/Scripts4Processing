clear 
close all, clc

% 
% Simple-West
load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/simple_west_M734.mat');
% Graphical Representation
f=figure('PaperType','a4')
f.Units='centimeters';
f.Position = [0.1 0.1 23 29]; 

% Spacing is always the same, calculate Lat, Lon only once
[X,Y] = xyz2LatLon(X,Y);

ax1=subplot(1,1,1);
sshafunction(X,Y,BT,WL,WH,2);  % 10 s --> 2nd element (now first to find Z indicies larger zero
title("t = 10s");

colorbarfkt()

% Load NetCDF file
ncFile = "/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Github/GMT/Iceland/input/GEBCO_08_Dec_2021_074bdff0ecf9/gebco_2021_n68.0_s62.0_w-25.0_e-12.0.nc";
lat = ncread(ncFile, 'lat');
lon = ncread(ncFile, 'lon');
bathymetryData = ncread(ncFile, 'elevation');
bathymetryData = transpose(bathymetryData);
[LON, LAT] = meshgrid(lon, lat);

hold on
%contour(LON,LAT,bathymetryData,[0 0], 'k', 'LineWidth', 1.5)
contour(LON,LAT,bathymetryData,[0 0], '-k', 'LineWidth', 1)


% Plot syn gauges
x_syn = [-17.9964, -18.8688, -18.6089, -17.3681, -18.5063, -18.1218];
y_syn = [66.518, 66.182, 66.091, 66.039, 65.965, 65.722];
hold on
plot(x_syn,y_syn,'r.','MarkerSize',14)

% Plot syn gauges numbers
x_syn_nr = x_syn + [0.03, -0.09, -0.09, 0.03, -0.09, -0.09]
y_syn_nr = y_syn
nr = {'6', '1', '2', '5', '3', '4'}
hold on
text(x_syn_nr, y_syn_nr, nr, 'Color','red','FontSize',14)

% Plot Hypo Simple-East
data = readtable("/import/freenas-m-05-seissol/kutschera/HIWI/./SeisSol/Nucleation_locations/epicentre_simple_east.csv")
x_epi = data{:,1}; 
y_epi = data{:,2};
[X_epi,Y_epi] = xyz2LatLon(x_epi,y_epi);
hold on
plot(X_epi,Y_epi,'MarkerFaceColor','yellow','Marker','pentagram','MarkerSize',16,'MarkerEdgeColor','black')







% % Give common xlabel, ylabel and title to your figure
% han=axes(f,'visible','off'); 
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% ylabel(han,'Latitude');
% xlabel(han,'Longitude');

%exportgraphics(f,'ssha.pdf')