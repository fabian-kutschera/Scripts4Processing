% -*- coding: utf-8 -*-
function sshafunction(X,Y,BT,WH,WL,ind)

%Z=WL(:,ind)-abs(BT(:,ind))+WH(:,ind);

%BT(BT>0)=nan; %makes topo invisible
Z=BT(:,ind).*max(BT(:,ind),0)+WH(:,ind);

% %reshape
% x = reshape(X,[],1);
% y = reshape(Y,[],1);
% 
% %utmZone
% zone='27 W';
% A = repmat(zone,length(x),1);
% 
% %Conversion
% [Lat, Lon]=utm2deg(x,y,A); 
% 
% %reshape pt2
% Lat = reshape(Lat, sqrt(length(Lat)), []);
% Lon = reshape(Lon, sqrt(length(Lon)), []);
Z=reshape(Z,sqrt(length(Z)),[]);
% 
% %plot
% pcolor(Lon,Lat,Z)

pcolor(X,Y,Z)
shading flat

%Crameri
load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/CrameriColourmaps/ScientificColourMaps7/cork/cork.mat')
%colormap(cork)
ax1 = gca;
set(ax1,'Colormap',cork)
caxis([-0.6, 0.6])


xlim([-20 -16.5])
ylim([65.4 66.9])

yticks([65.5, 66, 66.5]);
yticklabels({'65.5°N','66°N','66.5°N'});

xticks([-19,-17]);
xticklabels({'19°W','17°W'});


% time = (ind-1)*10;
% subtitle(['t=',num2str(time),'s'])

ax=gca;
ax.Layer='top';
ax.TickLength=ax.TickLength*1.5;
ax.FontSize=12;
ax.TickDir='both' %'in'

end
