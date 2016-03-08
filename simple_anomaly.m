CASE=1
if CASE == 1
    % soumya
    filepath = '/data/flow2/soumya/mean_13.8/'; id_start = 12026; nfiles = 150
    filepattern = 'mean_gmm_%d_%d.raw'
    output_filepattern_mean = 'anomaly/mean_b%d_%d.raw'
    id_step = 25
    W=151;  % Note!! Original range
    H=71
    D=56;
    res = 5;
else CASE==2
    % original
    filepath = '/home/chenchu/volumes/TEMP/14.20_16.00_14.20x2_150523/';      id_start = 49401
%     filepath = '/home/chenchu/volumes/TEMP/13.80_141219/';     id_start = 6201
%     filepath = '/home/chenchu/volumes/TEMP/14.20_16.00/';    id_start = 20601
    filepattern = 's35_noinj.r2b%d.p3d.g%d_Pressure.raw'
    output_filepattern_mean = 'anomaly/anomaly_b%d_%d.raw'
    id_step = 25
    nfiles = 576
    W=151;  % Note!! Original range
    H=71
    D=56;
    res = 1;
end

nblocks = 36;
w1 = floor(W/res);
h1 = floor(H/res);
d1 = floor(D/res);
ncells = w1*h1*d1
% Z_TH = 2.8% 3.3296;
for t=1:nfiles
    block_data = zeros(ncells, nblocks);
    for b=1:nblocks
        if CASE==1
            filename = sprintf(strcat(filepath, filepattern), b-1, id_start+id_step*(t-1) )
        else
            filename = sprintf(strcat(filepath, filepattern), b, id_start+id_step*(t-1) )
        end
        fp = fopen(filename, 'rb');
        data = fread(fp, ncells, 'float32');    
        fclose(fp);
        block_data(:,b) = data;
    end
    
    %% anomaly analysis
    anomaly = zeros(ncells,nblocks);
    for i=1:ncells
        data = block_data(i,:);
        s = std(data);
        m = mean(data);
        if s==0
            anomaly(i,:) = zeros(1, nblocks);
            continue;
        end
        anomaly(i,:) = (data - m)/s;            
    end

    for b=1:nblocks
        filename = strcat(filepath, sprintf(output_filepattern_mean, b, t) )
        fp = fopen(filename, 'wb')
        fwrite(fp, anomaly(:,b), 'float32');
        fclose(fp);
    end
    
end

