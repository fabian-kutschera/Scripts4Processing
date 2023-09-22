function plot_epi(dir)
data = readtable(dir);
x_epi = data{:,1};
y_epi = data{:,2};
[X_epi,Y_epi] = xyz2LatLon(x_epi,y_epi);
plot(X_epi,Y_epi,'MarkerFaceColor','yellow','Marker','pentagram','MarkerSize',15,'MarkerEdgeColor','black')

end