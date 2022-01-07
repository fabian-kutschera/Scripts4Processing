clear 
close all, clc


%% Graphical Representation
f=figure('PaperType','a4')
f.Units='centimeters';
%f.Position = [100 100 540 400];
f.Position = [0.1 0.1 23 29]; 

% Simple West
load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/simple_west_M734.mat');

% Spacing is always the same, calculate Lat, Long only once
[X,Y] = xyz2LatLon(X,Y);

ax1=subplot(6,3,1); % 10 s --> 2nd element
sshafunction(X,Y,BT,WL,WH,2);
title("t = 10s");
set(gca,'Xticklabel',[]);
ax2=subplot(6,3,2); % 120 s --> 13th element
sshafunction(X,Y,BT,WL,WH,13);
title("t = 120s");
set(gca,'Yticklabel',[]); 
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
%title('Simple West');
ax3=subplot(6,3,3); % 600 --> 61st element
sshafunction(X,Y,BT,WL,WH,61);
title("t = 600s");
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
colorbarfkt()

% Simple Middle
%load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/simple_middle_M733.mat');
ax4=subplot(6,3,4); % 10 s --> 2nd element
sshafunction(X,Y,BT,WL,WH,2);
set(gca,'Xticklabel',[]);
ax5=subplot(6,3,5); % 120 s --> 13th element
sshafunction(X,Y,BT,WL,WH,13);
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
%title('Simple Middle');
ax6=subplot(6,3,6); % 600 --> 61st element
sshafunction(X,Y,BT,WL,WH,61);
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
colorbarfkt()

% Simple East
%load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/simple_east_M734.mat');
ax7=subplot(6,3,7); % 10 s --> 2nd element
sshafunction(X,Y,BT,WL,WH,2);
set(gca,'Xticklabel',[]);
ax8=subplot(6,3,8); % 120 s --> 13th element
sshafunction(X,Y,BT,WL,WH,13);
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
%title('Simple East');
ax9=subplot(6,3,9); % 600 --> 61st element
sshafunction(X,Y,BT,WL,WH,61);
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
colorbarfkt()

% Complex West
%load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/complex_west_M674.mat');
ax10=subplot(6,3,10); % 10 s --> 2nd element
sshafunction(X,Y,BT,WL,WH,2);
set(gca,'Xticklabel',[]);
ax11=subplot(6,3,11); % 120 s --> 13th element
sshafunction(X,Y,BT,WL,WH,13);
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
%title('Complex West');
ax12=subplot(6,3,12); % 600 --> 61st element
sshafunction(X,Y,BT,WL,WH,61);
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
colorbarfkt()

% Complex Middle
%load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/complex_middle_M707.mat');
ax13=subplot(6,3,13); % 10 s --> 2nd element
sshafunction(X,Y,BT,WL,WH,2);
set(gca,'Xticklabel',[]);
ax14=subplot(6,3,14); % 120 s --> 13th element
sshafunction(X,Y,BT,WL,WH,13);
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
%title('Complex Middle');
ax15=subplot(6,3,15); % 600 --> 61st element
sshafunction(X,Y,BT,WL,WH,61);
set(gca,'Yticklabel',[]);
set(gca,'Xticklabel',[]); %to just get rid of the numbers but leave the ticks.
colorbarfkt()

% Complex East
%load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/complex_east_M668.mat');
ax16=subplot(6,3,16); % 10 s --> 2nd element
sshafunction(X,Y,BT,WL,WH,2);
ax17=subplot(6,3,17); % 120 s --> 13th element
sshafunction(X,Y,BT,WL,WH,13);
set(gca,'Yticklabel',[]); 
%title('Complex East');
ax18=subplot(6,3,18); % 600 s--> 61st element
sshafunction(X,Y,BT,WL,WH,61);
set(gca,'Yticklabel',[]); 
colorbarfkt()


% % Give common xlabel, ylabel and title to your figure
% han=axes(f,'visible','off'); 
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% ylabel(han,'Latitude');
% xlabel(han,'Longitude');

%exportgraphics(f,'ssha.pdf')