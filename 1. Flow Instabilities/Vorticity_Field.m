function [] = Vorticity_Field()
%VORTICITY_FIELD 
    strvar = 'linear';
    map=redblue(256);
    for i = 20:10:110
        clearvars -except i map strvar
        disp(i)
        load(['1030_',num2str(i),'.mat'])
        
        x = [ans.X];
        y = [ans.Y];
        u = [ans.U];%horizontal vel.
        v = [ans.V];%vertical vel.
        
        %{
        %1010-specific parameters (chosen based on velocity fields):
        xi = linspace(125,600,476); %x-range over which we interpolate the velocity-values to a regular grid 
        yi = linspace(20,565,546); %y-range over which we interpolate the velocity-values to a regular grid 
        %}
        %{
        %1020:
        xi = linspace(100,630,531);
        yi = linspace(20,565,546);  
        %}
        %%{
        %1030:
        xi = linspace(100,600,501);
        yi = linspace(20,565,546);  
        %}

        %Interpolating the velocity components over regular grids
        [x_grid,y_grid] =  meshgrid(xi,yi);
        u_gd = griddata(x,y,u,x_grid,y_grid,strvar);
        v_gd = griddata(x,y,v,x_grid,y_grid,strvar);
        
        %The spatial gradients of the velocity components:
        %gradient() only works as long as the spacing in each direction is 1, 
        % otherwise would have to input spacing as an argument
        [dvdx, dvdy] = gradient(v_gd); 
        [dudx, dudy] = gradient(u_gd);
        vort = dvdx - dudy; %vorticity is the curl of the velocity field

        %Plotting the vorticity (without smoothing)
        %{
        figure;
        colormap(map)
        camlight right;
        surf(x_grid,y_grid,vort,'EdgeColor','None');
        lighting phong;
        lim=caxis*0.1;
        lim=(abs(lim(1))+abs(lim(2)))/2;
        caxis([-lim, lim]);
        title('Vorticity plot');
        xlabel('x');
        ylabel('y');
        view(2)
        colorbar
        saveas(gcf,['Vort_',num2str(i),'_1030.png'])
        close all
        %}
        
        %Smoothing with convolution to smooth 2D data that contains high frequency components:
        for kvar = [3,5,10,50]
            K = (1/(kvar^2))*ones(kvar);
            Zsmooth1 = conv2(vort,K,'same');
                %Plotting
            figure;
            colormap("redblue")
            camlight right;
            surf(x_grid,y_grid,Zsmooth1,'EdgeColor','None')
            lighting phong;
            %keep colour axis centered on zero (white)
            lim=caxis;
            lim=(abs(lim(1))+abs(lim(2)))/2;
            caxis([-lim, lim]);
            title(['Vorticity plot (smoothed by convolving with a ',num2str(kvar),'x',num2str(kvar),' kernel)'])
            xlabel('x')
            ylabel('y')
            view(2)
            saveas(gcf,['Vort_',num2str(i),'_1030_Kern',num2str(kvar),'.png'])
            close all
        end
    end
end
