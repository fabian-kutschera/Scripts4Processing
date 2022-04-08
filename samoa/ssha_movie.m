clear all, close all, clc

load('/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Script4Processing/samoa/simple_middle_M733.mat');
%%
[X,Y] = xyz2LatLon(X,Y);
%%
f=figure;
f.Units='centimeters';
f.Position=[18 15 18 13];
set(f,'DefaultAxesFontSize',14)

movieid=VideoWriter('simple_middle_M733/ssha_simple_middle_M733');
movieid.FrameRate=4;
movieid.Quality=100;
open(movieid)


for ind=linspace(1,241,241)
    time = (ind-1)*10;
    sshafunction(X,Y,BT,WL,WH,ind);
    subtitle(sprintf('Tsunami simple Middle @ t = %4d s',time),'FontSize',14);
    xlabel('Longitude','FontSize',13);
    ylabel('Latitude','FontSize',13);
    
    h=colorbar();
    colorTitleHandle = get(h,'Title');
    titleString = 'ssha [m]';
    set(colorTitleHandle ,'String',titleString);
    h.Ticks = [-0.6 -0.4 -0.2 0.0 0.2 0.4 0.6];

%   print('movie/ssha.jpeg','-djpeg','-r600')
    print(sprintf('simple_middle_M733/ssha_t%d.jpeg', time),'-djpeg','-r600')


    %writeVideo(movieid,getframe(f));
    filename=sprintf('simple_middle_M733/ssha_t%d.jpeg', time);
    writeVideo(movieid,imread(filename));  
    %writeVideo(movieid,imread('movie/ssha.jpeg'));  

end

close(movieid)
