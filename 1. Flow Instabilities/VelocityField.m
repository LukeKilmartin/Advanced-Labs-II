clear all
figure()
tic
PredictiveTracker_open("C:\Users\Luke\Desktop\Advanced Laboratories I & II\Adv Labs 2\1 Flow Instabilities\tiff files\1010\20\*.tiff", ...
    10 ... %threshold brightness
    ,1 ... %within maxdisp particles of the predicted location
    ,         "C:\Users\Luke\Desktop\Advanced Laboratories I & II\Adv Labs 2\1 Flow Instabilities\tiff files\1010\Averages\20avg.tiff" ...
    ,1 ... %minarea size of particles
    ,0) %invert==0; PredictiveTracker seeks particles brighter than the background




%{
%structs:
%If any value input is an empty cell array, {}, then the output is an empty structure array. 
% %To specify a single empty field, use [].
x=[ans.X];
%x1={ans.X};
y=[ans.Y];
u=[ans.U]; %horizontal vel.
v=[ans.V]; %vertical vel.

%quiver(X,Y,U,V)
quiver(x,y,u,v)
toc
%}


