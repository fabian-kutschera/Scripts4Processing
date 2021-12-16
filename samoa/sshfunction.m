function sshfunction(X,Y,BT,WH,WL,ind)

Z=WL(:,ind)-abs(BT(:,ind))+WH(:,ind);
Z=reshape(Z,size(X,1),size(X,2));
pcolor(X,Y,Z)
shading flat

colormap(jet); % continuous colormap
%colormap(jet(10)) % discrete colormap with 10 values

caxis([-0.5, 0.5])

%colorbar

time = (ind-1)*10;
title(['t=',num2str(time),'s'])

ax=gca;
ax.Layer='top';
ax.TickLength=ax.TickLength*2;

end