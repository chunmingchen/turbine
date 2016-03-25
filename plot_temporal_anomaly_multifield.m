function gui_blending(varin)
    SAVE_IMG = 0
    global degmat1g degmat2g
%     global DISP_Y_RES SHOW_CASE  stepi filepath file_pattern TH res nfiles period W H D Ws Hs Ds endi degmat starti startb
    
    DISP_Y_RES = 6
    SHOW_CASE = 'PASSAGE'
    stepi = 1
    res = 1;

%     nfiles = 100 
    nfiles = 1438
    nfiles = 1436 % temporal
    period = 360
    W=30
    H=14
    D=11
    d1 = ceil(D/DISP_Y_RES);
    Ws = 1:W;
    Hs = 1:H;
    Ds = 1:d1;
    stepi = 10;
    startb = 1;

    RECALL_LABEL = 1

    
%     label = 'spatial pressure 14.2all'; filepath = {'/TEMP/uncertain/14.2x1/spatial_pressure_anomaly/out/','/TEMP/uncertain/14.2x2/spatial_pressure_anomaly/out/'}; file_pattern = 'emd_b%d_%d.raw'; starti = {20612, 35012}; TH = 0.45;
%     label = 'temporal pressure 14.2all'; filepath = {'/TEMP/uncertain/14.2x1/temporal_anomaly_pressure_raw/','/TEMP/uncertain/14.2x2/temporal_anomaly_pressure_raw/'}; file_pattern = 'temp_anomaly_b%d_%d.raw'; starti = {1,1}; stepi=1; TH = 0.15; startb=0;
    label = 'spatial pressure 16.0'; filepath = '/TEMP/uncertain/16.0/spatial_pressure_anomaly/raw/'; file_pattern = 'emd_b%d_%d.raw'; starti = 35022; TH = 0.45; nfiles = 717;
%     label = 'temporal pressure 16.0'; filepath = '/TEMP/uncertain/16.0/temporal_anomaly_pressure_raw/'; file_pattern = 'temp_anomaly_b%d_%d.raw'; starti = 1; stepi = 1; TH = .15; nfiles = 717; startb=0;
%--old
%     label = 'spatial pressure 14.2x1'; filepath = '/TEMP/uncertain/14.2x1/spatial_pressure_anomaly/out/'; file_pattern = 'emd_b%d_%d.raw'; starti = 20612; TH = 0.45;
%     label = 'spatial pressure 14.2x2'; filepath ='/TEMP/uncertain/14.2x2/spatial_pressure_anomaly/out/'; file_pattern = 'emd_b%d_%d.raw'; starti = 35012; TH = 0.45;

    load_data()

    %----------------------------------------------
    degmat1 = degmat;
    degmat1g = degmat;

    %----------------------------------------------

    RECALL_LABEL = 0
    
%     label = 'spatial entropy 14.2all'; filepath = {'/TEMP/uncertain/14.2x1/spatial_entropy_anomaly/raw/','/TEMP/uncertain/14.2x2/spatial_entropy_anomaly/raw/'};  file_pattern = 'emd_b%d_%d.raw'; starti = {20612, 35012}; TH = 0.35;
%     label = 'temporal entropy 14.2all'; filepath = {'/TEMP/uncertain/14.2x1/temporal_anomaly_entropy_raw/','/TEMP/uncertain/14.2x2/temporal_anomaly_entropy_raw/'}; file_pattern = 'temp_anomaly_b%d_%d.raw'; starti = {1, 1}; stepi = 1; TH = 0.12; startb=0;
    label = 'spatial entropy 16.0'; filepath = '/TEMP/uncertain/16.0/spatial_entropy_anomaly/raw/'; file_pattern = 'emd_b%d_%d.raw'; starti = 35022; TH = 0.35; nfiles = 717;
%     label = 'temporal entropy 16.0'; filepath = '/TEMP/uncertain/16.0/temporal_anomaly_entropy_raw/'; file_pattern = 'temp_anomaly_b%d_%d.raw'; starti = 1; stepi = 1; TH = 0.06; nfiles = 717; startb=0;
%--old
%     label = 'spatial entropy 14.2x1'; filepath = '/TEMP/uncertain/14.2x1/spatial_entropy_anomaly/raw/';  file_pattern = 'emd_b%d_%d.raw'; starti = 20612; TH = 0.35;
%     label = 'spatial entropy 14.2x2'; filepath = '/TEMP/uncertain/14.2x2/spatial_entropy_anomaly/raw/';  file_pattern = 'emd_b%d_%d.raw'; starti = 35012;  TH = 0.35;

    load_data()

    degmat2 = degmat;
    degmat2g = degmat;
    
% fout=fopen('test.raw','wt')
% fwrite(fout, degmat, 'float32');
% fclose(fout)

    %----------------------------------------------
    %% start plot
    fig = figure('Visible','off', 'Position', [20 800 1000 600], 'Color', [1 1 1], ...
        'SizeChangedFcn', @sizeChanged);

    txt = uicontrol('Style','text',...
        'Position',[15 15  200 20],'BackgroundColor','white', ...
        'String','Blending:', 'FontSize', 24);
    
    % slider
    sld = uicontrol('Style', 'slider',...
        'Min',0,'Max',1,'Value',0.5);
    addlistener(sld,'ContinuousValueChange', @surfzlim);
    fig.UserData = sld;
    sld.UserData = fig;
    
    sizeChanged(fig, 0);
        
    first_draw = true;
    surfzlim(0,0) % draw
    first_draw = false;
    
    set(fig,'Visible','on')
    
    %----------------------------------------------
    %% functions
    function sizeChanged(fig, callbackdata)
        pos = get(fig, 'Position');
        w = pos(3); h = pos(4);
        set(fig.UserData, 'Position', [w*.2 15 w*.75 20]); % x0 y0 w h
%         set(gca, 'Position', [20 0 1000 600])
        set(gca, 'Position', [.09 .21 .9 .77 ])
        set(txt, 'Position', [15 5  200 40]);
    end
   
    function surfzlim(source,callbackdata)
        set(fig, 'Visible', 'off');
        if ~first_draw
            L = get(gca,{'xlim','ylim'});  % Get axes limits.
        end
        val = get(sld,'value');        
        
        draw(val)
%         title('Anomaly Chart')

        if ~first_draw
            zoom reset
            set(gca,{'xlim','ylim'},L)
        end
        set(fig, 'Visible', 'on');
    end

    function cmap = tune_cmap(cmap)        
        cmap = cmap([1 5:64],:);
        cmap = [ [1 1 1] ; cmap ];
        %cmap = flipud(cmap);
    end

    function draw(alpha)
        degmat1(1,1) = 1000;
        degmat2(1,1) = 1000;
        maxval1 = min(3, max(max(degmat1)));
        maxval2 = min(3, max(max(degmat2)));
%         img1 = ind2rgb(ceil((degmat1/maxval1).^.9*70), tune_cmap(othercolor('Reds4')) );
%         img2 = ind2rgb(ceil((degmat2/maxval2).^.9*70), tune_cmap(othercolor('Bu_10')) );
        coexist = (degmat1/maxval1).*(degmat2/maxval2)*10*alpha*(1-alpha); % spatial
%         coexist = (degmat1/maxval1).*(degmat2/maxval2)*60*alpha*(1-alpha); % temporal
        img1 = ind2rgb(ceil((degmat1/maxval1 + coexist).^.8*70), tune_cmap(othercolor('Bu_10')) );
        img2 = ind2rgb(ceil((degmat2/maxval2 + coexist).^.8*70), tune_cmap(othercolor('Reds4')) );
        
        img = img1.*((1-alpha)) + img2.*(alpha);
        h = image(img);
        
%         alpha_data1 = (degmat1 > 0);
%         alpha_data2 = (degmat2 > 0);
%         h = image(img1);
%         set(h, 'AlphaData', alpha_data1*(1-alpha).^.3)
%         hold on
% %         
%         h = image(img2);
%         set(h, 'AlphaData', alpha_data2*(alpha).^.3)

        % xrange = starti+period-1:endi+period-1; % temporal anomaly
        xrange = 1:nfiles;
        yrange = 1:length(Ds)*36;

        maxval = max(max(degmat));
        minval = min(min(degmat));
        seg = (maxval-minval)/45;
        %h = imagesc(xrange, yrange, degmat, [minval-seg, (maxval)]);
        yticklabel = 1:2:36;
        %xticklabel = starti:stepi:endi;
        % set(gca,'YTickLabel',{'9' '8' '7' '6' '5' '4' '3' '2' '1' '0'})

        set(gca,'YDir','normal',  ...
            'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
            'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
            'ytick', 1:length(Ds)*2:(36*length(Ds)-1), ...
            'yticklabel', yticklabel, ...
            ... %'xticklabel', xticklabel, ...
            ... %'xtick', 1:skip:timesteps, ...
            'color', 'w', ...
            'YDir', 'reverse', ...
            'XColor', [.1 .1 .1],...
            'YColor', [.1 .1 .1] ,   ...
            'FontSize', 24, ...
            'GridLineStyle','-' ...
            )
        xlabel('Time Step')
        ylabel('Passage')

        hold off
        
        if 0
            out_image = sprintf('%s/anomaly.eps', '.') %filepath)
            saveas(gca, out_image, 'psc2')
        end
    end

    function load_data()
        mat_filename = strcat(label, '.mat');
        if RECALL_LABEL
            try
                load(mat_filename);
                disp(strcat(mat_filename, ' loaded'))
                return
            catch
                disp(strcat(mat_filename, ' not found'))
            end
        end
        
        if iscell(filepath)            
            for i = 1:length(filepath)                
                degmat_single = load_data_single(filepath{i}, starti{i});
                if i==1
                    degmat = degmat_single;
                else
%                     degmat = [degmat degmat(:,end-1:end) degmat_single];
                    degmat = [degmat degmat(:,end-3:end) degmat_single];
                    %                ^^^^^^^^^^^^^^^ skipped time steps
                end
            end
        else
            degmat = load_data_single(filepath, starti);
        end
        
        save(mat_filename, 'degmat');        
        disp(strcat(mat_filename, ' saved'))
    end
    
    function degmat = load_data_single(filepath, starti)
        
        endi = starti+stepi*(nfiles-1);
        DEGS = length(Ds)*36;
        degmat = zeros(DEGS, length(starti:stepi:endi));
        for i=1:nfiles
            if mod(i,100)==1
                disp(i)
            end
            for b=1:36
                filename = sprintf(strcat(filepath, file_pattern), b+startb-1, starti+stepi*(i-1))
                fp = fopen(filename, 'rb');
                data = fread(fp, W*H*D, 'float32');
                fclose(fp);

                list = find(abs(data)>TH);
                for j=1:length(list)
                    idx = list(j)-1;
                    xx = mod(idx, length(Ws))+1;
                    yy = mod(floor(idx/length(Ws)), length(Hs))+1;            
                    zz = floor(idx/length(Ws)/length(Hs))+1;
                    zz = floor((zz-1)/DISP_Y_RES)+1; %!!!!!!!!!!
                    y = zz+ length(Ds)*(b-1);

                    degmat(y, i) = degmat(y, i) + 1;
                end
            end

        end

    end
end


