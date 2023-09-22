clear all, close all, clc

% Load NetCDF file once
ncFile = "/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Github/GMT/Iceland/input/GEBCO_08_Dec_2021_074bdff0ecf9/gebco_2021_n68.0_s62.0_w-25.0_e-12.0.nc";
lat = ncread(ncFile, 'lat');
lon = ncread(ncFile, 'lon');
bathymetryData = ncread(ncFile, 'elevation');
bathymetryData = transpose(bathymetryData);
[LON, LAT] = meshgrid(lon, lat);

f = figure();
%[ha, pos] = tight_subplot(6,3,[.01 .03],[.04 .04],[.08 .1])
[ha, pos] = tight_subplot(6,3,[.01 .03],[.04 .04],[.11 .11])
%[ha, pos] = tight_subplot(rows,columns,[vspace_between_plots hspace_between_plots],[bottom_margin top_margin],[left_margin right_margin]); 
set(gcf,'Units','centimeters');
%set(gcf,'Position',[0.1 0.1 19.5 26]);
set(gcf,'Position',[0.1 0.1 18.5 22]);

ind1 = 2;
time1 = (ind1-1)*10;
ind2 = 13;
time2 = (ind2-1)*10;
ind3 = 61;
time3 = (ind3-1)*10;

% Simple West
load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/simple_west_M734.mat');
% Spacing is always the same, calculate Lat, Long only once
[X,Y] = xyz2LatLon(X,Y);
axes(ha(1)); 
sshafunction(X,Y,BT,WL,WH,ind1);
subtitle(['t=',num2str(time1),'s']);
text(-20,66.8,'Simple-West','Units','data','Color','[0.6350, 0.0780, 0.1840]','FontSize',12,'FontWeight','bold');
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)
% Plot syn gauges
x_syn = [-17.9964, -18.8688, -18.6089, -17.3681, -18.5063, -18.1218];
y_syn = [66.518, 66.182, 66.091, 66.039, 65.965, 65.722];
hold on
plot(x_syn,y_syn,'r.','MarkerSize',12)
hold on
% Plot syn gauges numbers
x_syn_nr = x_syn + [0.12, -0.31, -0.31, 0.12, -0.31, -0.31]
y_syn_nr = y_syn
nr = {'6', '1', '2', '5', '3', '4'}
text(x_syn_nr, y_syn_nr, nr, 'Color','red','FontSize',10)
hold on
% Plot Hypo Simple-West
plot_epi("/import/freenas-m-05-seissol/kutschera/HIWI/SeisSol/Nucleation_locations/epicentre_simple_west.csv")

axes(ha(2)); 
sshafunction(X,Y,BT,WL,WH,ind2);
subtitle(['t=',num2str(time2),'s']);
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)
axes(ha(3)); 
sshafunction(X,Y,BT,WL,WH,ind3);
subtitle(['t=',num2str(time3),'s']);
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)

% Simple Middle
load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/simple_middle_M733.mat');
% Spacing is always the same, calculate Lat, Long only once
[X,Y] = xyz2LatLon(X,Y);
axes(ha(4)); 
sshafunction(X,Y,BT,WL,WH,ind1);
text(-20,66.8,'Simple-Middle','Units','data','Color','[0.6350, 0.0780, 0.1840]','FontSize',12,'FontWeight','bold');
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)
hold on
% Plot Hypo Simple-Middle
plot_epi("/import/freenas-m-05-seissol/kutschera/HIWI/SeisSol/Nucleation_locations/epicentre_simple_middle.csv")

axes(ha(5)); 
sshafunction(X,Y,BT,WL,WH,ind2);
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)
axes(ha(6)); 
sshafunction(X,Y,BT,WL,WH,ind3);
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)

% Simple East
load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/simple_east_M734.mat');
% Spacing is always the same, calculate Lat, Long only once
[X,Y] = xyz2LatLon(X,Y);
axes(ha(7)); 
sshafunction(X,Y,BT,WL,WH,ind1);
text(-20,66.8,'Simple-East','Units','data','Color','[0.6350, 0.0780, 0.1840]','FontSize',12,'FontWeight','bold');
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)
hold on
% Plot Hypo Simple-East
plot_epi("/import/freenas-m-05-seissol/kutschera/HIWI/SeisSol/Nucleation_locations/epicentre_simple_east.csv")

axes(ha(8)); 
sshafunction(X,Y,BT,WL,WH,ind2);
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)
axes(ha(9)); 
sshafunction(X,Y,BT,WL,WH,ind3);
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)

% Complex West
load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/complex_west_M674.mat');
% Spacing is always the same, calculate Lat, Long only once
[X,Y] = xyz2LatLon(X,Y);
axes(ha(10)); 
sshafunction(X,Y,BT,WL,WH,ind1);
text(-20,66.8,'Complex-West','Units','data','Color','[0.6350, 0.0780, 0.1840]','FontSize',12,'FontWeight','bold');
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)
hold on
% Plot Hypo Complex-West
plot_epi("/import/freenas-m-05-seissol/kutschera/HIWI/SeisSol/Nucleation_locations/complex_locations/epicentre_complex_west.csv")

axes(ha(11)); 
sshafunction(X,Y,BT,WL,WH,ind2);
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)
axes(ha(12)); 
sshafunction(X,Y,BT,WL,WH,ind3);
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)

% Complex Middle
load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/complex_middle_M707.mat');
% Spacing is always the same, calculate Lat, Long only once
[X,Y] = xyz2LatLon(X,Y);
axes(ha(13)); 
sshafunction(X,Y,BT,WL,WH,ind1);
text(-20,66.8,'Complex-Middle','Units','data','Color','[0.6350, 0.0780, 0.1840]','FontSize',12,'FontWeight','bold');
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)
hold on
% Plot Hypo Complex-Middle
plot_epi("/import/freenas-m-05-seissol/kutschera/HIWI/SeisSol/Nucleation_locations/complex_locations/epicentre_complex_middle.csv")

axes(ha(14)); 
sshafunction(X,Y,BT,WL,WH,ind2);
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)
axes(ha(15)); 
sshafunction(X,Y,BT,WL,WH,ind3);
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)

% Complex East
load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/complex_east_M668.mat');
% Spacing is always the same, calculate Lat, Long only once
[X,Y] = xyz2LatLon(X,Y);
axes(ha(16)); 
sshafunction(X,Y,BT,WL,WH,ind1);
text(-20,66.8,'Complex-East','Units','data','Color','[0.6350, 0.0780, 0.1840]','FontSize',12,'FontWeight','bold');
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)
hold on
% Plot Hypo Complex-East
plot_epi("/import/freenas-m-05-seissol/kutschera/HIWI/SeisSol/Nucleation_locations/complex_locations/epicentre_complex_east.csv")

axes(ha(17)); 
sshafunction(X,Y,BT,WL,WH,ind2);
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)
axes(ha(18)); 
sshafunction(X,Y,BT,WL,WH,ind3);
hold on
contour(LON,LAT,bathymetryData,[0 0], 'k-', 'LineWidth', 0.4)

set(ha(1:15),'XTickLabel',''); 
set(ha(2:3),'YTickLabel','')
set(ha(5:6),'YTickLabel','')
set(ha(8:9),'YTickLabel','')
set(ha(11:12),'YTickLabel','')
set(ha(14:15),'YTickLabel','')
set(ha(17:18),'YTickLabel','')

h=colorbar();
pos=h.Position;
hx=pos(1)*1.1; hy=pos(2)*8.9; hw=pos(3)*0.8; hh=pos(4)*2;
%hx=pos(1)*1.08; hy=pos(2)*8.9; hw=pos(3)*0.8; hh=pos(4)*2;
h.Position=[hx, hy, hw, hh];
h.Ticks = [-0.6 -0.4 -0.2 0.0 0.2 0.4 0.6];

ylabel(h, 'ssha [m]', 'Rotation', 0);
h.Label.Position(2)=h.Label.Position(2)+0.75;
h.Label.Position(1)=h.Label.Position(1)-3; %x-direction

grid off

exportgraphics(f,'fig07.pdf','Resolution', 300)
print('fig07','-dpng','-r450')