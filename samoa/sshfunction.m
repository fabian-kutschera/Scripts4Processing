function sshfunction(X,Y,BT,WH,WL,ind)

%Z=WL(:,ind)-abs(BT(:,ind))+WH(:,ind);
Z=BT(:,ind).*max(BT(:,ind),0)+WH(:,ind);

%Z=reshape(Z,size(X,1),size(X,2));
Z=reshape(Z,sqrt(length(Z)),[]);

pcolor(X,Y,Z)
shading flat

%colormap(jet); % continuous colormap
%colormap(jet(10)) % discrete colormap with 10 values
%colorbar
load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/CrameriColourmaps/ScientificColourMaps7/cork/cork.mat')
colormap(cork)
%colorbar

caxis([-0.5, 0.5])


time = (ind-1)*10;
subtitle(['t=',num2str(time),'s'])

ax=gca;
ax.Layer='top';
ax.TickLength=ax.TickLength*2;
ax.FontSize=14;
ax.TickDir='in' %'both'

end