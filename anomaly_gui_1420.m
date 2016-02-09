function gui(varin)
    global gdata
    SHOW_CASE = 'PASSAGE'
%     SHOW_CASE = 'ANGLE' % 0% user speed
%     SHOW_CASE = 'USER'
    USER_SPEED_PERCENT =  45
    SHOW_HOUGH = 0
    SHOW_IMLINE = 1
    CASE=1
    SAVE_IMG=1
    xextract = 7:10;
    fontsize=30; % 12
    timesteps1 = 576
    target_path1 = '/data/flow2/turbine_Stg/s35_noinj_14.20_150112_turb_6201-20601/output_important/anomaly/'
    timesteps2 = 576
    target_path2 = '/data/flow2/turbine_Stg/s35_noinj_14.20_150112_turb_continue_150216_20601-35001/saved/anomaly'
    timesteps3 = 576
    target_path3 = '/data/flow2/turbine_Stg/s35_noinj_14.20_150112_turb_continue_150216_continue_150301_35001-/anomaly'
    
    if nargin==0
        var = 'Pressure'
%         var = 'Entropy'
%         var = 'TotalPressure'
%         var = 'Density'
%         var = 'Temperature'
%         var = 'VelocityMagnitude'
    else
        var = varin
    end
    
    data_p1 = myload_data(var, timesteps1, target_path1);
    data_p2 = myload_data(var, timesteps2, target_path2);
    data_p3 = myload_data(var, timesteps3, target_path3);
    timesteps = timesteps1+timesteps2+timesteps3;
    data_p = struct;
    data_p.v = cat(3, data_p1.v, data_p2.v, data_p3.v)
%     data_p.v = zeros(size(data_p1,1), size(data_p1,2), timesteps);
%     data_p.v(:,:,1:timesteps1) = data_p1.v;
%     data_p.v(:,:,(timesteps1+1):timesteps2) = data_p2;
%     data_p.v(:,:,(timesteps2+1):timesteps3) = data_p3;
    
    gdata = data_p

    
    
    % Create a figure and axes
    f = figure('Visible','off', 'Position', [20 70 1800 1000]);

    % colormap(hot)
    cmap = colormap(hot);
    cmap = cmap([1:44,64],:);
    cmap = flipud(cmap);
    colormap(cmap)
    
    % Create pop-up menu
%     popup = uicontrol('Style', 'popup',...
%            'String', {'parula','jet','hsv','hot','cool','gray'},...
%            'Position', [20 340 100 50],...
%            'Callback', @setmap);    
    
   % Create push button
%     btn = uicontrol('Style', 'pushbutton', 'String', 'Clear',...
%         'Position', [20 20 50 20],...
%         'Callback', 'cla');       

    % Add a text uicontrol to label the slider.
    
   % Create slider
    if ~SAVE_IMG
        sld = uicontrol('Style', 'slider',...
            'Min',0,'Max',1000,'Value',0,...
            'Position', [220 45 1400 20],...
            'Callback', @surfzlim); 
        
        txt = uicontrol('Style','text',...
            'Position',[0 45 160 20],...
            'String','Saliency in size:');    
    end

    % Add a text uicontrol to label the slider.
%     txt = uicontrol('Style','text',...
%         'Position',[0 20 160 20],...
%         'String','Rotation speed percentage:');

    % checkbox for hough    
    if ~SAVE_IMG
        chk = uicontrol('Style','checkbox','String','Line Detection', ...
                       'Value',SHOW_HOUGH,'Position',[0 20 160 20],        ...
                        'Callback',@surfzlim);
    end
    
    % Make figure visble after adding all components
    surfzlim(0,0)
    set(f,'Visible','on')
    

    function setmap(source,callbackdata)
        val = source;
        maps = source.String;
        newmap = maps{val};
        colormap(newmap);
    end

    function surfzlim(source,callbackdata)
        if ~SAVE_IMG
            SHOW_HOUGH = get(chk, 'value');
            val = get(sld,'value');
        else 
            val = 0;
        end
        
       %subplot(1,2,1)
       show(data_p, (val/1000).^2)
       title(strcat(var, ' Anomaly'))
       %subplot(1,2,2)
       %show(data_d, val/1000)
       %title('Density Anomaly')
    end
    
    function data = myload_data(var, timesteps, target_path)
        data1=zeros(36*56, 16, timesteps);
        for t=1:timesteps
            file=sprintf('%s/Salient_%s_iso3.3296_%d.txt', target_path, var, t-1);
%             file=sprintf('%s/Salient_%s_iso3.5_%d.txt', target_path, var, t-1);
            data=dlmread(file);
            
            if strcmp(SHOW_CASE, 'ANGLE') || strcmp(SHOW_CASE, 'USER')
                if strcmp(SHOW_CASE, 'USER')
                    ratio = 1-USER_SPEED_PERCENT/100;
                else
                    ratio = 1;
                end
                passage_len = length(data)/36;
                rotate = mod( round((t-1)*passage_len/4*ratio) , length(data));
                if rotate>0
                    data = [data((rotate+1):end, :); data(1:rotate, :)];
                end
            end
            
            data1(:,:,t)=data;
        end
        data = struct;
        data.v = data1;
         
        data.v = data.v(:, xextract, :);
        
% 
%         % Determine the connected components:
%             conn = 6
%         data.CC = bwconncomp(data.v, conn);
%         data.S = regionprops(data.CC, 'Area');
%         data.S1 = regionprops(data.CC, data.v, 'MeanIntensity');
%         % Remove small objects:
%         data.L = labelmatrix(data.CC);
%         data.maxsum = max([data.S.Area].*[data.S1.MeanIntensity])
    end

    function show(sdata, ratio)
        if 0
            % remove small area
            P = 500
            P1 = sdata.maxsum*ratio 
            % P1= 3406207*.5  %entropy
            % P1 = 19885237 *.95  %velocity
            % P1 = 273366 *.05 %density
            % Compute the area of each component:
            BW2 = ismember(sdata.L, find(([sdata.S.Area].*[sdata.S1.MeanIntensity]) >= P1));
            data1 = sdata.v .*BW2;
        else
            data1 = sdata.v;
        end
        render(data1)

        if  SHOW_HOUGH
           show_hough(data1)
        end

        if SHOW_IMLINE
            h = imline(gca, [10 100], [100 100])
            id = addNewPositionCallback(h,@showspeed);
        end
                
        if SAVE_IMG
%             filename = strcat(target_path, '/anomaly.eps')

%           set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 6])
            filename = 'anomaly.eps'
            saveas(gca, filename, 'psc2')
        end
    end


    function render(data1)        
        degrees = size(data1,1);
        disp(size(data1,2))
        data2 = squeeze(sum(data1,2));
        
        I = find(data2>0);
        
        minval = min(min(min(data2(I))))
        
        maxval = max(max(max(data2)))
        
        seg = (maxval-minval)/44;
%         imagesc(sqrt(data2), [0 sqrt(maxval)])
%         imagesc(data2>0, [0 40])
        imagesc((data2), [minval-seg (maxval)])
        
        if SAVE_IMG
            skip = 200 %50
        else
            skip = 10
        end
        set(gca,'YDir','normal',  ...
            'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
            'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
            'ytick', 1:56:(36*56-1), ...
            ... 'xtick', 1:skip:timesteps, ...
            'color', 'w', ...
            'YDir', 'reverse', ...
            'XColor', [.5 .5 .5],...
            'YColor', [.5 .5 .5] ,   ...
            'FontSize', fontsize ...
            )
        if strcmp(SHOW_CASE, 'ANGLE')
            set(gca, 'yticklabel', 0:10:360);
            ylabel('Degree')
            
        elseif strcmp(SHOW_CASE, 'PASSAGE')
            set(gca, 'yticklabel', 1:36);
            ylabel('Passage')
            
        elseif strcmp(SHOW_CASE, 'USER')
            set(gca, 'yticklabel', 1:36);
            ylabel('')
        end
        
        ylim([1, degrees])
        xlim([1, timesteps])
        
        hold on
        
        for i=1:4
            if strcmp(SHOW_CASE, 'PASSAGE')
                x=[(i-1)*144+1; i*144+1];
                y=[1; 36*56+1];
                plot([x(1);x(2)], [y(1);y(2)], ':', 'color', [.5 .5 .5])
            elseif strcmp(SHOW_CASE, 'ANGLE')
                x=[(i-1)*144+1; i*144+1];
                y=[1; 36*56+1];
                plot([x(2);x(1)], [y(1);y(2)], '--', 'color', [.5 .5 .5])
            elseif strcmp(SHOW_CASE, 'USER')
                rate = 144/(1-USER_SPEED_PERCENT/100);
                x=[(i-1)*rate+1; i*rate+1];
                y=[1; 36*56+1];
                plot([x(2);x(1)], [y(1);y(2)], '--', 'color', [.5 .5 .5])
            end
        end
        
        hold off
        
        cbr = colorbar;
        xlabel('Time Step')
        set(cbr, 'YTick', [0])
        
    end

    function show_hough(data)
        BW = squeeze(sum(data, 2));
%         BW = [BW;BW]; % make repletive
%         rotI = imrotate(I,33,'crop');
%         BW = edge(rotI,'canny');
%         [H,T,R] = hough(BW);
        [H, T, R] = weightedHough(BW);
%         figure
%         imagesc(H);
%         pause
%         imshow(H,[],'XData',T,'YData',R,...
%                     'InitialMagnification','fit');
%         imagesc(H)
        
        axis on, axis normal, hold on;
        peaks = 30000;  ratio = 0.5;
        P  = houghpeaks(H,peaks,'threshold',ceil(ratio*max(H(:))));
        x = T(P(:,2)); y = R(P(:,1));
        plot(x,y,'s','color','white');
        % Find lines and plot them
        lines = houghlines(BW,T,R,P,'FillGap',20,'MinLength',10);
        hold on
        max_len = 0;
        
        roter_spd = 360/(3600/25);
        %disp(sprintf('Rotor speed: %f degrees/timestep', roter_spd))
        
        for k = 1:length(lines)
           xy = [lines(k).point1; lines(k).point2];
           plot(xy(:,1),xy(:,2),'--','LineWidth',1,'Color','black');
%            plot(xy(:,1),xy(:,2)-size(BW,1),'--','LineWidth',2,'Color','green');

           % Plot beginnings and ends of lines
%            plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%            plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

           % Determine the endpoints of the longest line segment
           len = norm(lines(k).point1 - lines(k).point2);
           if ( len > max_len)
              max_len = len;
              xy_long = xy;
           end
           spd = roter_spd - ((xy(2,2)-xy(1,2))/(56*36)*360)/(xy(2,1)-xy(1,1));
           percent = spd / roter_spd *100;
           disp(sprintf('Rotation speed: %.1f rotor speed, Length: %f', percent, len))
           str1 = sprintf('%.1f%%', percent);
           text(mean(xy(:,1)),mean(xy(:,2)),str1)
        end

        % highlight the longest line segment
        if length(lines)>0
            plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','blue');
%             plot(xy_long(:,1),xy_long(:,2)-size(BW,1),'LineWidth',2,'Color','blue');
        end
        hold off
    end

    function [h, T, R] = weightedHough(img)
        n = floor(2*norm(size(img)));
        RhoResolution = 1;
        D = sqrt((n - 1)^2 + (180 - 1)^2);
        diagonal = ceil(D/2/RhoResolution); 
        nrho = 2*diagonal + 1;
        T = -60:.5:-5; %-90:89
        ntheta = length(T);
        h = zeros(nrho, ntheta);
        for y=1:size(img,1)
            for x=1:size(img,2)
                v = img(y,x);
                if v>0
                    for th = T 
                        rho = round(x*cos(th/180*pi)+y*sin(th/180*pi));
                        h(rho+diagonal+1, th*2+121) = h(rho+diagonal+1, th*2+121) + v;
                    end
                end
            end
        end
        R = -diagonal:diagonal;
    end

    texthandle = 0;
    
    function showspeed(pos)
        roter_spd = 360/(3600/25);
        % pos is a 2-by-2 array [X1 Y1; X2 Y2].
        spd = roter_spd - ((pos(2,2)-pos(1,2))/(56*36)*360)/(pos(2,1)-pos(1,1)); 
        percent = spd / roter_spd *100;
        disp(sprintf('Rotation speed: %.1f rotor speed', percent))
        str1 = sprintf('%.1f%%', percent);
        if texthandle~=0
            try
                delete (texthandle)
            catch
            end
        end
        texthandle = text(mean(pos(:,1)),mean(pos(:,2)),str1);
    end
end