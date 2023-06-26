%Load the right file
clearvars
i=20;
load(['1010_',num2str(i),'.mat'])

x = [ans.X];
y = [ans.Y];
u = [ans.U]; %horizontal vel.
v = [ans.V]; %vertical vel.
%curl(x,y,u,v) - won't work - need to griddata x,y,u,v first
%min(x)
%max(x)
%min(y)
%max(y)

%manually determined for each image set (1010/1020/1030):
%%{
%1010-specific parameters (chosen based on velocity fields):
xi = linspace(125,600,476); %x-range over which we interpolate the velocity-values to a regular grid 
yi = linspace(20,565,546); %y-range over which we interpolate the velocity-values to a regular grid 
%}
%{
%1020:
xi = linspace(100,630,531);
yi = linspace(20,565,546);  
%}
%{
%1030:
xi = linspace(100,600,501);
yi = linspace(20,565,546);  
%}

strvar = 'linear'

[x_grid,y_grid] = meshgrid(xi,yi);
u_gd = griddata(x,y,u,x_grid,y_grid,strvar);
v_gd = griddata(x,y,v,x_grid,y_grid,strvar);

%gradient() only works as long as the spacing in each direction is 1, 
% otherwise would have to input spacing as an argument
[dvdx, dvdy] = gradient(v_gd); 
[dudx, dudy] = gradient(u_gd);

map=redblue(256);

vort = dvdx - dudy;

%%{
%Plotting the vorticity (without smoothing)
figure;
colormap(map)
camlight right;
surf(x_grid,y_grid,vort,'EdgeColor','None');
lighting phong;
lim=caxis*0.1;
lim=(abs(lim(1))+abs(lim(2)))/2;
caxis([-lim, lim]);
title(['Vorticity at 20 mA']);
xlabel('x');
ylabel('y');
zlabel('ω');
colorbar
view(2) %view in 2D from above - see only x,y axes
%}

%{
%Plotting the velocity interpolations onto a regular grid:
    %Horizontal velocity:
figure
mesh(x_grid,y_grid,u_gd)
colormap(jet)
camlight right
lighting phong
title(['U (Interpolation method == ',strvar,')'])
xlabel('x')
ylabel('y')
zlabel('Horizontal velocity (U)')

    %Vertical Velocity
figure
mesh(x_grid,y_grid,v_gd)
colormap("jet")
camlight right
lighting phong
title(['V (Interpolation method == ',strvar,')'])
xlabel('x')
ylabel('y')
zlabel('Vertical velocity (V)')
%}

%{
%Angular velocity using a curl black box:
[curlz,c_ang_vel]=curl(x_grid,y_grid,u_gd,v_gd);
    %Plotting:
figure;
colormap(map);
camlight right;
surf(x_grid,y_grid,c_ang_vel,'EdgeColor','None')
lighting phong;
caxis manual;
caxis([-.2 .2]);
title(['ang\_vel\_gd on peaks surface (',strvar,')'])
xlabel('x')
ylabel('y')
view(2)
%}

%%{
%Smoothing with convolution:
for kvar = [50];
    K = (1/(kvar))*ones(kvar);
    Zsmooth1 = conv2(vort,K,'same');
        %Plotting
    figure;
    colormap("redblue")
    camlight right;
    surf(x_grid,y_grid,Zsmooth1,'EdgeColor','None')
    lighting phong;
    lim=caxis;
    lim=(abs(lim(1))+abs(lim(2)))/2;
    caxis([-lim, lim]);
    title(['SMOOTH (',num2str(kvar),') vort on surface (',strvar,')'])
    xlabel('x')
    ylabel('y')
    colorbar
    view(2)
end
%}

%{
%The curlz:
figure
colormap("jet")
camlight right
lighting phong
surf(x_grid,y_grid,curlz,'EdgeColor','None')
title 'curlz on peaks surface'
xlabel('x')
ylabel('y')
view(2)
%}