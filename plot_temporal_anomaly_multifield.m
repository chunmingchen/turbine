function gui_blending(varin)
    SAVE_IMG = 0

%     global DISP_Y_RES SHOW_CASE  stepi filepath file_pattern TH res nfiles period W H D Ws Hs Ds endi degmat starti startb
    DISP_Y_RES = 6
    SHOW_CASE = 'PASSAGE'
    stepi = 1
    res = 1;

%     nfiles = 300 
    nfiles = 1437 %1339
    period = 360
    W=30
    H=14
    D=11
    d1 = ceil(D/DISP_Y_RES);
    Ws = 1:W;
    Hs = 1:H;
    Ds = 1:d1;
    stepi = 10
    startb = 1


    
    %filepath =
    %'/home/chenchu/volumes/TEMP/uncertain/gmm_files_generated_pressure/anomaly/';     starti=6212    
%     filepath = '/home/chenchu/volumes/TEMP/uncertain/14.2x1/pressure_anomaly/out/';  starti = 21612
    filepath = '/home/chenchu/volumes/TEMP/uncertain/14.2x2/pressure_anomaly/out/';  starti = 35012
%     filepath = '/home/chenchu/volumes/TEMP/uncertain/14.2x1/entropy_anomaly/raw/';  starti = 21612
%     filepath = '/home/chenchu/volumes/TEMP/uncertain/14.2x2/entropy_anomaly/raw/';  starti = 35012
    file_pattern = 'emd_b%d_%d.raw'
    TH = 0.45 %3.3296

    
    load_data()

    %----------------------------------------------
    degmat1 = degmat;

    %----------------------------------------------

    
    %filepath =
    %'/home/chenchu/volumes/TEMP/uncertain/gmm_files_generated_pressure/anomaly/';     starti=6212    
    %     filepath = '/home/chenchu/volumes/TEMP/uncertain/14.2x1/pressure_anomaly/out/';  starti = 21612
    %     filepath = '/home/chenchu/volumes/TEMP/uncertain/14.2x2/pressure_anomaly/out/';  starti = 35012
    %     filepath = '/home/chenchu/volumes/TEMP/uncertain/14.2x1/entropy_anomaly/raw/';  starti = 21612
%     filepath = '/home/chenchu/volumes/TEMP/uncertain/14.2x2/entropy_anomaly/raw/';  starti = 35012
%     TH = 0.35 %3.3296

    % temporal:
    filepath = '/TEMP/uncertain/temporal_anomaly/14.2_x2_raw/temporal_anomaly_pressure/'; starti = 0; stepi = 1;
%     filepath = '/TEMP/uncertain/temporal_anomaly/14.2_x1_raw/temporal_anomaly_pressure/'; starti = 0; stepi = 1;
    file_pattern = 'temp_anomaly_b%d_%d.raw';
    TH = 0.22
    startb = 0

    load_data()

    degmat2 = degmat;

% fout=fopen('test.raw','wt')
% fwrite(fout, degmat, 'float32');
% fclose(fout)

    %----------------------------------------------
    %% plot
    f = figure('Visible','off', 'Position', [20 70 1800 1000], 'Color', [1 1 1]);
    
    % slider
    sld = uicontrol('Style', 'slider',...
        'Min',0,'Max',1,'Value',0.5,...
        'Position', [220 45 1400 20],...
        'Callback', @surfzlim); 

    
    surfzlim(0,0)
    set(f,'Visible','on')
   
    function surfzlim(source,callbackdata)
        if ~SAVE_IMG
            val = get(sld,'value');
        else 
            val = 0;
        end
        disp(sprintf('alpha: %f', val) )
        
        draw(val)
        title('Anomaly Chart')

    end

    function cmap = tune_cmap(cmap)        
        cmap = cmap([1:64],:);
        %cmap = flipud(cmap);
    end

    function draw(alpha)
        
        maxval = max(max(degmat1));
        img1 = ind2rgb(ceil((degmat1/maxval).^.7*100), tune_cmap(othercolor('Reds4')) );
        alpha_data1 = (degmat1 > 0);
        maxval = max(max(degmat2));
        img2 = ind2rgb(ceil((degmat2/maxval).^.7*100), tune_cmap(othercolor('Bu_10')) );
        alpha_data2 = (degmat2 > 0);
        
        h = image(img1);
        set(h, 'AlphaData', alpha_data1*(1-alpha).^.3)
        hold on
        
        h = image(img2);
        set(h, 'AlphaData', alpha_data2*(alpha).^.3)

        % xrange = starti+period-1:endi+period-1; % temporal anomaly
        xrange = 1:nfiles;
        yrange = 1:length(Ds)*36;

        maxval = max(max(degmat));
        minval = min(min(degmat));
        seg = (maxval-minval)/45;
        %h = imagesc(xrange, yrange, degmat, [minval-seg, (maxval)]);
        yticklabel = 1:36;
        %xticklabel = starti:stepi:endi;
        % set(gca,'YTickLabel',{'9' '8' '7' '6' '5' '4' '3' '2' '1' '0'})

        set(gca,'YDir','normal',  ...
            'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
            'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
            'ytick', 1:length(Ds):(36*length(Ds)-1), ...
            'yticklabel', yticklabel, ...
            ... %'xticklabel', xticklabel, ...
            ... %'xtick', 1:skip:timesteps, ...
            'color', 'w', ...
            'YDir', 'reverse', ...
            'XColor', [.3 .3 .3],...
            'YColor', [.3 .3 .3] ,   ...
            'FontSize', 12, ...
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
        endi = starti+stepi*(nfiles-1)
        DEGS = length(Ds)*36;
        degmat = zeros(DEGS, length(starti:stepi:endi));
        for i=1:nfiles
            if mod(i,100)==1
                disp(i)
            end
            for b=1:36
                filename = sprintf(strcat(filepath, file_pattern), b+startb-1, starti+stepi*(i-1));
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


