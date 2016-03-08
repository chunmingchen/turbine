function gui(varin)
    global gdata
    SHOW_CASE = 'PASSAGE'
%     SHOW_CASE = 'ANGLE' % 0% user speed
%     SHOW_CASE = 'USER'
    USER_SPEED_PERCENT =  45
    SHOW_HOUGH = 0
    SHOW_IMLINE = 1
    CASE=16
    SAVE_IMG=1
%     xextract = 1:16 
    xextract = 7:10;
    COLORING = 0 % count
%     COLORING = 1 % depth
    switch (CASE)
        case 16
            timesteps = 576
            target_path = '~/volumes/TEMP/s35_noinj_14.20_16.00_14.20x2_150523/output'
                    
    end 
    if nargin==0
%         var = {'Entropy'};
%         var = {'Pressure'};
%         var = {'TotalPressure'};
%         var = {'Density'};
%         var = {'VelocityMagnitude'};
        var = {'Entropy','Pressure','TotalPressure','Density','VelocityMagnitude'};
            %,'Temperature',
    else
        var = varin
    end
    
    data_p = cell(length(var),1);
    for i=1:length(var)
        data_p{i} = myload_data(var{i});
    end
%     data_d = myload_data('Density');
    
    gdata = data_p

    
    
    % Create a figure and axes
    f = figure('Visible','off', 'Position', [20 70 1800 1000], 'Color', [1 1 1]);

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
    set(f,'Visible','on');
    

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
        
        show(data_p, (val/1000).^2)
%         title(strcat(var{1}, ' Anomaly'))
        %show(data_d, val/1000)
%title('Density Anomaly')
    end
    
    function data = myload_data(var)

        if iscell(target_path) % multiple files
            ts = 0;
            for i=1:length(timesteps)
                ts = ts + timesteps{i};
            end
            data1=zeros(36*56, 16, ts);
            ts = 1;
            for i=1:length(timesteps)
                disp(sprintf('loading %s...', target_path{i}));
                data1(:,:,ts:ts+timesteps{i}-1) = myload_data1(var, timesteps{i}, target_path{i});
                ts = ts + timesteps{i};
            end
            timesteps = ts-1;
        else
            data1 = myload_data1(var, timesteps, target_path);
        end
        data = struct;
        data.v = data1(:, xextract, :);        

        % Determine the connected components:
            conn = 6
        data.CC = bwconncomp(data.v, conn);
        data.S = regionprops(data.CC, 'Area');
        data.S1 = regionprops(data.CC, data.v, 'MeanIntensity');
        % Remove small objects:
        data.L = labelmatrix(data.CC);
        data.maxsum = max([data.S.Area].*[data.S1.MeanIntensity])
    end

    function data1 = myload_data1(var, timesteps, target_path)
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
    end

    function cmap = create_colormap(a,b)
        n = 44;
        R = linspace(a(1), b(1), n);
        G = linspace(a(2), b(2), n);
        B = linspace(a(3), b(3), n);
        R = [1 R];
        G = [1 G];
        B = [1 B];
        cmap = [R' G' B'];
    end

    texthandle = 0;
    function show(sdata, ratio)
        
        for i=1:length(sdata)
            n=64;
            if i==1  
                cmap = create_colormap([1 .7 .7], [1 0 0]);
            elseif i==2
                cmap = create_colormap([.9 .9 .3], [.7 .7 0]);
            elseif i==3
                cmap = create_colormap([.7 1 .7], [0 1 0]);
            elseif i==4
                cmap = create_colormap([.6 .9 .9], [0 .6 .6]);
            elseif i==5
                cmap = create_colormap([.7 .7 1], [0 0 1]);
            end
            colormap( cmap );
            
            render(sdata{i}.v, var{i})
            
            
            hold on
        end
        hold off

    end


    function render(data1, varname)        
        disp(size(data1,2))
        if COLORING==0
            data2 = squeeze(sum(data1,2));

            I = find(data2>0);
            minval = min(min(min(data2(I))))
            maxval = max(max(max(data2)))
        else
            data2 = zeros(size(data1,1), timesteps)+max(xextract)+1;
            for t=1:timesteps
                for i=1:size(data1,1)
                    idx = find(data1(i,:,t)>0);
                    if length(idx)>0
                        data2(i,t) = idx(1);
                    end
                end
            end
            minval = min(xextract);
            maxval = max(xextract);
        end

        
        seg = (maxval-minval)/44;
%         imagesc(sqrt(data2), [0 sqrt(maxval)])
%         imagesc(data2>0, [0 40])

        imagesc((data2), [minval-seg (maxval)])
        
%         g = data2 / maxval;
%         set(h, 'AlphaData', g);
        
        if strcmp(SHOW_CASE, 'ANGLE')
            yticklabel = 0:10:360;
            ylabel('Degree')
            
        elseif strcmp(SHOW_CASE, 'PASSAGE')
            yticklabel = 1:36;
            ylabel('Passage')
            
        elseif strcmp(SHOW_CASE, 'USER')
            yticklabel = 1:36;
            ylabel('')
        end
        degrees = size(data2,1);
        set(gca,'YDir','normal',  ...
            'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
            'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
            'ytick', 1:56:(36*56-1), ...
            'yticklabel', yticklabel, ...
            ... %'xtick', 1:skip:timesteps, ...
            'color', 'w', ...
            'YDir', 'reverse', ...
            'XColor', [.5 .5 .5],...
            'YColor', [.5 .5 .5] ,   ...
            'FontSize', 12 ...
            )
        ylim([1, degrees])
        xlim([1, timesteps])
        
        
        if SAVE_IMG
            skip = 50
        else
            skip = 10
        end
        
        hold on
        
        for i=1:(size(data2,2)/144)
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
        
        %cbr = colorbar;
        xlabel('Time Step')
        %set(cbr, 'YTick', [0])
        
        if  SHOW_HOUGH
           show_hough(data1)
        end

        if SHOW_IMLINE
            texthandle = imline(gca, [10 100], [100 100])
            id = addNewPositionCallback(texthandle,@showspeed);
        end
                
        if SAVE_IMG
%             filename = strcat(target_path, '/anomaly.eps')
            filename = strcat(varname, '_anomaly.eps')
            saveas(gca, filename, 'psc2')
        end
    end

%%---------------------------------------------------------------------------------------------
    %% OK functions
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

    %% Hough
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

end