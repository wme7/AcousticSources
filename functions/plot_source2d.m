function plot_source2d(varargin)
% Produce a visualization of the field in the (x,y)-plane

switch nargin
    case 6
        % Expected Inputs
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        field = varargin{4};
        name = varargin{5};
        range = varargin{6};

        % Get field size
        [~,~,nz] = size(field); 
        
        % Get x3-direction middle plane
        zc = round(nz/2);

        % Plot the field surface
        surf(x(:,:,zc),y(:,:,zc),field(:,:,zc));
        axis equal; view(2); colorbar; shading interp; 
        title(name,'Interpreter','latex');
        xlabel('$x_1$','Interpreter','latex'); 
        ylabel('$x_2$','Interpreter','latex');
        if not(isempty(range))
            caxis(range);
        else
            caxis([-1.0,1.0]);
        end

    case 11
        % Expected Inputs
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        field = varargin{4};
        x_s = varargin{5};
        y_s = varargin{6};
        z_s = varargin{7};
        T_s = varargin{8};
        field_s = varargin{9};
        name = varargin{10};
        range = varargin{11};

        % Get field size
        [~,~,nz] = size(field); 
        
        % Get x3-direction middle plane
        zc = round(nz/2);
        
        % Field
        subplot(221)
            surf(x(:,:,zc),y(:,:,zc),field(:,:,zc)); hold on
            scatter3(x_s,y_s,z_s,'.k'); hold off
            axis equal; view(2); colorbar; shading interp; 
            title(name,'Interpreter','latex');
            xlabel('$x_1$','Interpreter','latex'); 
            ylabel('$x_2$','Interpreter','latex');
            if not(isempty(range))
                caxis(range);
            else
                caxis([-1.0,1.0]);
            end
        subplot(223)
            imagesc(field_s);
            title('Registered Signal','Interpreter','latex'); 
            xlabel('capteurs','Interpreter','latex'); 
            ylabel('time [ms]','Interpreter','latex');
        subplot(2,2,[2,4])
            polarplot(T_s,rms(field_s)); thetatickformat('degrees')

    otherwise
        disp('Usage:')
        disp('plot_source2d(x,y,z,field,name)')
        disp('plot_source2d(x,y,z,field,xs,ys,zs,Ts,s_field,name)')
end

end % function