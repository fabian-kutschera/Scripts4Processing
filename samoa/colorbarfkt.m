function colorbarfkt()

h=colorbar('Location','east');
pos=h.Position;
hx=pos(1)*1.06; hy=pos(2); hw=pos(3)*0.5; hh=pos(4)*2;
h.Position=[hx, hy, hw, hh];
%set(h,'YTickLabel',{'-0.5'; '-0.1'; '0.0'; '0.1'; '0.5'});
%h.Ticks = [-1.0 -0.5 0.0 0.5 1.0];
h.Ticks = [-0.6 -0.4 -0.2 0.0 0.2 0.4 0.6];

% h.Label.Position = [pos(1)/2 pos(2)+2];
% h.Label.Rotation = 0;
% h.Label.String = ('ssha [m]');
% yt = get(h, 'YTick');
% set(h,'XTickLabel',sprintf('%2.4|', yt));
% % h.Label.HorizontalAlignment = 'right';


colorTitleHandle = get(h,'Title');
titleString = 'ssha [m]';
set(colorTitleHandle ,'String',titleString);
h.Label
h.Label.Position = [pos(1) pos(2)*2];
% 
figure_handle = gcf;
cbar_handle = findobj(figure_handle,'tag','Colorbar')
set(cbar_handle, 'YAxisLocation','right')

end