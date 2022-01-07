clear 
close all, clc


%% Graphical Representation
f=figure()
f.Units='centimeters';
%f.Position = [100 100 540 400];
%f.Position = [0.1 0.1 20 29]; 

% Simple West
load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/simple_west_M734.mat');
ax1=subplot(1,1,1); % 10 s --> 2nd element
sshafunction(X,Y,BT,WL,WH,2);
%colorbarfkt()

%% Coord trafo
load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/simple_west_M734.mat');

%% reshape

re_x = reshape(X,[],1);
re_y = reshape(Y,[],1);
%%
%A = zeros(length(re_x),1);
zone='27 W';
A = repmat(zone,length(re_x),1);

%%
[Lat, Lon]=utm2deg(re_x,re_y,A); 

%% reshape pt2

lat = reshape(Lat, sqrt(length(Lat)), []);
lon = reshape(Lon, sqrt(length(Lon)), []);


%% plot
figure(10)
ax1=subplot(1,1,1); % 10 s --> 2nd element
sshfunction(lon,lat,BT,WL,WH,2);
colorbarfkt()

xlim([-20 -16.5])
ylim([65.4 66.9])