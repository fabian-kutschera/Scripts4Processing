% Load coastline data
load coastlines

% Create a new figure
figure;

% Plot coastline
plot(coastlon, coastlat, 'k');

% Set axis limits
xlim([-180 180]);
ylim([-90 90]);

% Set aspect ratio
daspect([1 1 1]);

% Set grid
grid on;

% Set labels
xlabel('Longitude');
ylabel('Latitude');

% Set title
title('Coastline Plot');
%%
% Load NetCDF file
ncFile = "/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Github/GMT/Iceland/input/GEBCO_08_Dec_2021_074bdff0ecf9/gebco_2021_n68.0_s62.0_w-25.0_e-12.0.nc"

lat = ncread(ncFile, 'lat');
lon = ncread(ncFile, 'lon');
bathymetryData = ncread(ncFile, 'elevation');
bathymetryData = transpose(bathymetryData);

[LON, LAT] = meshgrid(lon, lat);

figure;
%imagesc(lon,lat,bathymetryData)
contour(LON,LAT,bathymetryData,[0 0], 'k')
