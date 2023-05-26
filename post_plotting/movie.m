movieid=VideoWriter('output_combined/ssvv');
movieid.FrameRate=4;
movieid.Quality=100;
open(movieid)


for time=linspace(0.0,179.9,1800)
    filename = sprintf('./output_combined/trace_%06.2fs.png', time)
    writeVideo(movieid,imread(filename));   

end

close(movieid)
