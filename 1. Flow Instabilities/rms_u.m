clearvars
strvar = 'linear';
rmsu = zeros(1,10);
rmsvel = zeros(1,10);
    for i = 20:10:110
        clearvars -except i rmsu strvar rmsvel
        %disp(i)
        load(['1010_',num2str(i),'.mat'])
        
        x = [ans.X];
        y = [ans.Y];
        u = [ans.U];%horizontal vel.
        v = [ans.V];%vertical vel.
        
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

        %Interpolating the velocity components over regular grids
        [x_grid,y_grid] =  meshgrid(xi,yi);
        u_gd = griddata(x,y,u,x_grid,y_grid,strvar);
        v_gd = griddata(x,y,v,x_grid,y_grid,strvar);

        vsqr=v_gd.^2;
        usqr=u_gd.^2;
        %total velocity:
        velocity_magnitude = sqrt(vsqr+usqr);
        velocity_mag_sqr = velocity_magnitude.^2;
        rms_velocity = sqrt(nanmean(velocity_mag_sqr(:)));
        rmsvel(i/10 - 1) = rms_velocity;

        %horizontal
        mean_usqr = nanmean(usqr(:));
        sqrt_mean_usqr = sqrt(mean_usqr);
        rmsu(i/10 - 1) = sqrt_mean_usqr;

    end
    scatter([20,30,40,50,60,70,80,90,100,110],rmsvel,"filled")
    ylabel("Reynolds Number",'Interpreter','latex');
    xlabel('Driving current $I$ ($mA$)','Interpreter','latex');
    saveas(gcf,'re_vs_I_1010.png')

    scatter(rmsvel,rmsu,"filled")
    ylabel("RMS x-velocity: $\left\langle u_{x}^{2}\right\rangle^{1 / 2}$",'Interpreter','latex');
    xlabel('Reynolds number','Interpreter','latex');
    saveas(gcf,'ux_vs_Re_1010.png')
    %{
    scatter([20,30,40,50,60,70,80,90,100,110],rmsu,"filled")
    ylabel("RMS x-velocity: $\left\langle u_{x}^{2}\right\rangle^{1 / 2}$",'Interpreter','latex');
    xlabel('Driving current $I$ ($mA$)','Interpreter','latex');
    saveas(gcf,['x_vel_1010.png'])
    %}
    close all
