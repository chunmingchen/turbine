CASE=1 
% CASE=2 % single passage

if CASE==1  % full annulus
    filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_13.80_150127_regular_sampled/raw/'
    filepattern = '14.2-13.8-turbine_%d.vti_Pressure_slice38.raw'
    output_filepattern = 'anomaly/slice38_%d.raw'
    id_step = 1
    id_start = 1
    nfiles = 576
    period = 144
    W=509
    H=509
else % single passage
    slice = '20'; %'40'
    filepath = '/media/chenchu/TEMP/14.20_13.8-0slice/'
    % filepattern = strcat('s35_noinj.r2b20.p3d.g%d_Pressure_slice',slice,'.raw')
    filepattern = strcat('s35_noinj.r2b10.p3d.g%d_Pressure_slice',slice,'.raw')
    %output_filepattern = strcat('tempdiff_slice',slice,'_%d.raw')
    output_filepattern = strcat('anomaly_b10_slice',slice,'_%d.raw')
    id_step = 25
    id_start = 0
    nfiles = 534
    period = 144
    W=71
    H=56
end

blocks = 36

points_per_block = period/blocks;


window = period;
first = true;
data = zeros(W, H, window);
endi = (nfiles-window);
for i=1:endi
    if i==1
        for j=i:i+window-1
            filename = sprintf(strcat(filepath, filepattern), id_start+id_step*(j-1) )

            fp = fopen(filename, 'rb');
            data(:,:,j) = fread(fp, [W,H], 'float32');
            fclose(fp);
        end
    else
        data(:,:,1:end-1) = data(:,:,2:end);
        
        filename = sprintf(strcat(filepath, filepattern), id_start+id_step*(i+window-1) )
        fp = fopen(filename, 'rb');
        data(:,:,end) = fread(fp, [W,H], 'float32');
        fclose(fp);
    end

    data_z = zeros(size(data));
    for p=1:points_per_block
        pts_id = p:points_per_block:window;
        data_sym = data(:,:,pts_id) ;
        dim_=3;
        std_mat = std(data_sym, 1, dim_);
        std_mat(std_mat==0) = 1;
        data_z(:,:,pts_id) = abs((data_sym - repmat(mean(data_sym, dim_), [1,1,blocks]) ) ./ repmat(std_mat, [1,1,blocks]));
    end
    
    if i==1
        startj = 1;
    else
        startj = floor(window/2);
    end
    if i==endi
        endj = window;
    else
        endj = floor(window/2);
    end
    for j=startj:endj
        filename = sprintf(strcat(filepath, output_filepattern), i+j-1)
        fp = fopen(filename, 'wb');
        fwrite(fp, data_z(:,:,j), 'float32');
        fclose(fp);

                imagesc(squeeze(data_z(:,:,j)))
                colorbar
                drawnow
    end
end
        


