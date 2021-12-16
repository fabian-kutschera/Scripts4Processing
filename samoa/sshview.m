clear all
close all, clc

load('/import/freenas-m-05-seissol/bo/CHEESE/forFabian/west.mat');

%% Graphical Representation

f=figure;
f.Units='centimeters';
f.Position=[8 13 35 9]

ax1=subplot(1,3,1); % 10 s --> 2nd element
sshfunction(X,Y,BT,WL,WH,2);

ax2=subplot(1,3,2); % 60 s --> 7th element
sshfunction(X,Y,BT,WL,WH,7);

ax3=subplot(1,3,3); % 600 --> 61st element
sshfunction(X,Y,BT,WL,WH,61);

h=colorbar('Location','east')
%h.Position
h.Label.String='ssh [m]'
