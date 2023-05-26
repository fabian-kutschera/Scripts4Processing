clear all
close all
tic
datadir='/import/freenas-m-05-seissol/kutschera/HIWI/samoa/outputSamoa_EGU2021/outputSamoa_tanioka_complex_west_M6.74/flash_20210315_161521.470_0_xdmf.h5';
datadir1='/import/freenas-m-05-seissol/kutschera/HIWI/samoa/outputSamoa_EGU2021/outputSamoa_tanioka_complex_west_M6.74/flash_20210315_161521.470_1_xdmf.h5';
datadir2='/import/freenas-m-05-seissol/kutschera/HIWI/samoa/outputSamoa_EGU2021/outputSamoa_tanioka_complex_west_M6.74/flash_20210315_161521.470_2_xdmf.h5';

geometry = hdf5read(datadir,'/20/g');
waterheight=hdf5read(datadir,'/20/a/k');
bathymetry=hdf5read(datadir,'/20/a/b');
%% TESTING
ind = 20;
[X,Y]=meshgrid(min(geometry(1,:)):150:max(geometry(1,:)),min(geometry(2,:)):150:max(geometry(2,:)));
geometry = hdf5read(datadir,['/',num2str(ind-1),'/g']);
waterheight=hdf5read(datadir,['/',num2str(ind-1),'/a/k']);
grid_waterheight=griddata(geometry(1,:),geometry(2,:),waterheight(:),X(:),Y(:));
%%
scatter(geometry(1,:),geometry(2,:),[],waterheight(:),'filled')
toc
%%
%set new uniform girds for all time steps
tic
%regrid with 150 m spacing
[X,Y]=meshgrid(min(geometry(1,:)):150:max(geometry(1,:)),min(geometry(2,:)):150:max(geometry(2,:)));
T=[-1; (10:10:2400)']; %time steps
% load('/import/freenas-m-05-seissol/bo/CHEESE/forFabian/west.mat');
WH=zeros(length(X(:)),length(T));
WL=WH;
BT=WH;
for ind=1:100
geometry = hdf5read(datadir,['/',num2str(ind-1),'/g']);
waterheight=hdf5read(datadir,['/',num2str(ind-1),'/a/k']);
grid_waterheight=griddata(geometry(1,:),geometry(2,:),waterheight(:),X(:),Y(:));
waterlevel=hdf5read(datadir,['/',num2str(ind-1),'/a/h']);
grid_waterlevel=griddata(geometry(1,:),geometry(2,:),waterlevel(:),X(:),Y(:));
bathymetry=hdf5read(datadir,['/',num2str(ind-1),'/a/b']);
grid_bathymetry=griddata(geometry(1,:),geometry(2,:),bathymetry(:),X(:),Y(:));
%%to check the regrid output matches with the paraview view
% grid_waterheight=reshape(grid_waterheight,size(X,1),size(X,2));
% figure (10)
% imagesc(X(1,:),Y(:,1),grid_waterheight)
% set(gca,'YDir','normal')
WH(:,ind)=grid_waterheight;
WL(:,ind)=grid_waterlevel;
BT(:,ind)=grid_bathymetry;
clear geometry waterheight grid* waterlevel bathymetry
end
toc
tic
for ind=101:200
geometry = hdf5read(datadir1,['/',num2str(ind-1),'/g']);
waterheight=hdf5read(datadir1,['/',num2str(ind-1),'/a/k']);
grid_waterheight=griddata(geometry(1,:),geometry(2,:),waterheight(:),X(:),Y(:));
waterlevel=hdf5read(datadir1,['/',num2str(ind-1),'/a/h']);
grid_waterlevel=griddata(geometry(1,:),geometry(2,:),waterlevel(:),X(:),Y(:));
bathymetry=hdf5read(datadir1,['/',num2str(ind-1),'/a/b']);
grid_bathymetry=griddata(geometry(1,:),geometry(2,:),bathymetry(:),X(:),Y(:));
% grid_waterheight=reshape(grid_waterheight,size(X,1),size(X,2));
% figure (10)
% imagesc(X(1,:),Y(:,1),grid_waterheight)
% set(gca,'YDir','normal')
WH(:,ind)=grid_waterheight;
WL(:,ind)=grid_waterlevel;
BT(:,ind)=grid_bathymetry;
clear geometry waterheight grid* waterlevel bathymetry
end
toc
tic
for ind=201:241
geometry = hdf5read(datadir2,['/',num2str(ind-1),'/g']);
waterheight=hdf5read(datadir2,['/',num2str(ind-1),'/a/k']);
grid_waterheight=griddata(geometry(1,:),geometry(2,:),waterheight(:),X(:),Y(:));
waterlevel=hdf5read(datadir2,['/',num2str(ind-1),'/a/h']);
grid_waterlevel=griddata(geometry(1,:),geometry(2,:),waterlevel(:),X(:),Y(:));
bathymetry=hdf5read(datadir2,['/',num2str(ind-1),'/a/b']);
grid_bathymetry=griddata(geometry(1,:),geometry(2,:),bathymetry(:),X(:),Y(:));
% grid_waterheight=reshape(grid_waterheight,size(X,1),size(X,2));
% figure (10)
% imagesc(X(1,:),Y(:,1),grid_waterheight)
% set(gca,'YDir','normal')
WH(:,ind)=grid_waterheight;
WL(:,ind)=grid_waterlevel;
BT(:,ind)=grid_bathymetry;
clear geometry waterheight grid* waterlevel bathymetry
end

%save complex_west_M674.mat WH WL BT X Y -v7.3
toc
