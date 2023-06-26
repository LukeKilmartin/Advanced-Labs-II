function [] = Velocity_Field()
%VELOCITY_FIELD Summary of this function goes here
%   Detailed explanation goes here


for i = 20:10:110
    clearvars -except i
    disp(i)
    load(['1020_',num2str(i),'.mat'])
    %structs:
    %If any value input is an empty cell array, {}, then the output is an empty structure array. 
    % %To specify a single empty field, use [].
    x=[ans.X];
    %x1={ans.X};
    y=[ans.Y];
    u=[ans.U]; %horizontal vel.
    v=[ans.V]; %vertical vel.
    disp()
    %quiver(x,y,u,v) %vector plot
    %saveas(gcf,['V_F_',num2str(i),'.png'])
end
%clear all
close all
end