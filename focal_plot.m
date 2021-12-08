clear all
close all

addpath(genpath('./MTplot'))

datadir=('./output/');

%fnames=dir([datadir '*.h5']);
fnames=dir([datadir 'PointSourceFile_1_1_Simple_West_M7.343.h5']);
focal_info=zeros(length(fnames),6);
for istn=1:length(fnames)
    hinfo=hdf5info([datadir fnames(istn).name]);
Mt=hdf5read([datadir,fnames(istn).name],'/MomentTensor');
fm=Mt;
 M(1,1) = fm(2); M(2,2) = fm(3); M(3,3) = fm(1);
        M(2,1) = -fm(6); M(1,2) = M(2,1);
        M(3,1) = fm(4); M(1,3) = M(3,1);
        M(3,2) = -fm(5); M(2,3) = M(3,2);

[Strike1,Dip1,Rake1,Strike2,Dip2,Rake2]=MT2SDR(M)
focal_info(istn,:)=[Strike1,Dip1,Rake1,Strike2,Dip2,Rake2];
MTplot(M)
titlename=strcat('SDR1:',num2str(round(Strike1*100)/100),';',...
    num2str(round(Dip1*100)/100),';',num2str(round(Rake1*100)/100),';SDR2:',...
    num2str(round(Strike2*100)/100),';',num2str(round(Dip2*100)/100),';',num2str(round(Rake2*100)/100));
title(titlename,'Fontsize',18)
set(gca,'LineWidth',2)
set(gca,'Fontsize',15,'Fontweight','bold')
name0=strsplit(fnames(istn).name,'.')
print(gcf,'-djpeg',[datadir name0{1},'_focal'])
%saveas(gcf, [datadir 'Simple_West_M7.343.png'])

%clear Mt M Dip* Strike* Rake* titlename name0
%close all
end