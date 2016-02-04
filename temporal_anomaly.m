% CASE=1 
% CASE=1.5
% CASE=2 % single passage
CASE=3 % single passage all points
% CASE=3.1 % 14.20x2
RUN_TEMPDIFF = 0;
RUN_BOUNDS = 0;
RUN_PCA = 1;
RUN_LI2016 = 0;
if CASE==1  % full annulus
    filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_13.80_150127_regular_sampled/raw/'
    filepattern = '14.2-13.8-turbine_%d.vti_Pressure_slice38.raw'
    output_filepattern = 'tempdiff_slice38_%d.raw'
    id_step = 1
    id_start = 1
    nfiles = 576
    period = 144
    W=509
    H=509
    D=1
    Ws=1:W
    Hs=1:H
    Ds=1:D
elseif CASE==1.5 % full annulus     low res
    filepath = '/media/chenchu/TEMP/14.20_13.80_raw_low/'
    filepattern = 'sampling_res0.005_%d.vti_Pressure.raw'
    id_step = 25
    id_start = 0
    nfiles = 534
    period = 144
    W=35
    H=204
    D=204
    Ws=1:W
    Hs=1:H
    Ds=1:D
    
%     SLICE=16; % leading edge
%     output_filepattern = 'tempdiff/tempdiff_slice16_%d.raw'
%     output_filepattern_LI2016 = 'li2016/li2016_slice16_%d.raw'
    % no slicing
    output_filepattern = 'tempdiff/tempdiff_%d.raw'
    output_filepattern_LI2016 = 'li2016/li2016_%d.raw'
elseif CASE==2 % single passage
    slice = '20'; %'40'
    filepath = '/media/chenchu/TEMP/14.20_13.8_single/'
    % filepattern = strcat('s35_noinj.r2b20.p3d.g%d_Pressure_slice',slice,'.raw')
    filepattern = strcat('s35_noinj.r2b10.p3d.g%d_Pressure_slice',slice,'.raw')
    %output_filepattern = strcat('tempdiff_slice',slice,'_%d.raw')
    output_filepattern = strcat('tempdiff_b10_slice',slice,'_%d.raw')
    id_step = 25
    id_start = 0
    nfiles = 534
    period = 144
    W=71
    H=56
    D=1
    Ws=1:W
    Hs=1:H
    Ds=1:D
elseif CASE==3 % single passage from full annulus
%     block = 10
    filepath = '/media/chenchu/TEMP/14.20_13.80_single/'
    filepattern = strcat('s35_noinj.r2b',int2str(block),'.p3d.g%d_Pressure.raw')
    output_filepattern = strcat('tempdiff/tempdiff_b',int2str(block),'_%d.raw')
    id_start = 0
    id_step = 25
    nfiles = 534
    period = 36
    W=151;  % Note!! Original range
    H=71
    D=56
    Ws = 11:2:81;
    Hs = 40:2:H;
    Ds = 1:2:D;
elseif CASE==3.1 % single passage from full annulus
%     block = '10'
    filepath = '/home/chenchu/volumes/TEMP/14.20_16.00_14.20x2_150523/'
    filepattern = strcat('s35_noinj.r2b',int2str(block),'.p3d.g%d_Pressure.raw')
    output_filepattern = strcat('tempdiff/tempdiff_b',int2str(block),'_%d.raw')
    id_start = 49401
    id_step = 25
    nfiles = 576
    period = 36
    W=151;  % Note!! Original range
    H=71
    D=56
    Ws = 11:2:81;
    Hs = 40:2:H;
    Ds = 1:2:D;
end

%%
window = period;
if RUN_LI2016
    window = period*2;
end
test = [];
test_fft = [];
test_dist = [];
n=length(Ws)*length(Hs)*length(Ds);
data = zeros(n, window);
ts = nfiles-window;
freqs = floor(window/2)+1;
if RUN_BOUNDS
    freq_max=zeros(freqs,1);
    freq_mean=zeros(freqs,1);
end
for i=1:ts
    if i==1
        for j=i:i+window-1
            filename = sprintf(strcat(filepath, filepattern), id_start+id_step*(j-1) )
            fp = fopen(filename, 'rb');
            a = fread(fp, [W*H*D,1], 'float32');

            a = reshape(a, [W,H,D]);
            a = a(Ws,Hs,Ds);
            a = reshape(a, [n,1]);
            
            data(:,j) = a;
            fclose(fp);
        end
    else
        %data = circshift(data, -1, 2);
        data(:,1:end-1) = data(:,2:end);
        
        filename = sprintf(strcat(filepath, filepattern), id_start+id_step*(i+window-2) )
        fp = fopen(filename, 'rb');
        a = fread(fp, [W*H*D,1], 'float32');
        
        a = reshape(a, [W,H,D]);
        a = a(Ws,Hs,Ds);
        a = reshape(a, [n,1]);

        data(:,end) = a;
        fclose(fp);
    end
    %test(i) = data(88,400,window);

    if RUN_PCA || RUN_TEMPDIFF || RUN_BOUNDS
        mat = zeros(floor(window/2)+1, n);
        for x=1:n
            s = squeeze(data(x,:));

            %display(sprintf('%d %d', ((t-1)*W+1) ,(t*W)))
            Y = fft(s);
            P2 = abs(Y)/window; % Fourier amplitude
            P1 = P2(1:window/2+1);
            P1(2:end-1) = 2*P1(2:end-1);

            mat(:,x) = P1;
        end
    end
    
    %% RUN_TEMPDIFF
    if RUN_TEMPDIFF
        test_fft(:,i) = P1;
        if i~=1
            % distance
            diff = mat-premat;
            dist = sqrt(sum(diff.^2,1));

            filename = sprintf(strcat(filepath, output_filepattern), i);
            
            fp = fopen(filename, 'wb');
            fwrite(fp, dist, 'float32');
            fclose(fp);

            test_dist(i) = dist(1,x);
        end
        premat = mat;
    end % if RUN_TEMPDIFF
    
    %% RUN_BOUNDS
    if RUN_BOUNDS
        freq_max = max(freq_max, max(mat, [], 2));
        freq_mean = freq_mean + mean(mat, 2);
    end
    
    %% RUN_PCA
    if RUN_PCA
        if i==1
            nzidx = find(sum(mat,1)>0);
            nzcols = length(nzidx) % nonzero columns
            coeffs = zeros(freqs, nzcols*ts, 'single');
        end
        coeffs(:,((i-1)*nzcols+1):(i*nzcols)) = mat(:,nzidx);
    end
    
    %% RUN_LI2016
    if RUN_LI2016
        %mat = data(:, period+1:period*2-1); % prev period
        %premat = data(:, 1:period); % cur period
        range = 4;
        range1 = 1:range;
        range2 = period+1:period+range;
        
        matmag = sqrt(sum( data(:,range2).^2,2));
        prematmag = sqrt(sum(data(:,range1).^2,2));
        Rc = -dot(data(:,range1), data(:,range2), 2)./matmag./prematmag;

        filename = sprintf(strcat(filepath, output_filepattern_LI2016), i);

        fp = fopen(filename, 'wb');
        fwrite(fp, Rc, 'float32');
        fclose(fp);
        
    end
    
end

if RUN_BOUNDS
    freq_mean = freq_mean / ts;
    figure; plot(freq_max);
    figure; plot(freq_mean);
    save(strcat(filepath, 'freq_max.mat'), 'freq_max');
    save(strcat(filepath, 'freq_mean.mat'), 'freq_mean');
end
if RUN_PCA
    freq_mean = mean(coeffs, 2, 'native');
    coeffs = coeffs - repmat(freq_mean, 1,size(coeffs,2));
    
    F = 1/(nzcols-1) * coeffs*coeffs';
    save(strcat(filepath, 'F.mat'), 'F');
    save(strcat(filepath, 'nzidx.mat'), 'nzidx');
    [PC, V] = eig(F);
    V=diag(V);
    [junk, I] = sort(-1*V);
    V  = V(I);       % V: sorted large to small
    PC = PC(:,I);
    PC1 = PC(:, 1:5);
    projected = PC1' * coeffs;
%     save(strcat(filepath, 'coeffs.mat'), 'coeffs');
%     projected1 = reshape(projected, [2, nzcols, ts]);
    % recons = PC1 * projected;
    
    fp = fopen(strcat(filepath, 'projected.raw'), 'wb');
    fwrite(fp, projected, 'float32');
    fclose(fp);
    
    fp = fopen(strcat(filepath, 'projected.txt'), 'wt');
    fprintf(fp, '%d %d %d\n', length(Ws), length(Hs), length(Ds) );
    
    fprintf(fp, '%d %d %d\n', size(projected,1), nzcols, ts);  % h-d coefficients
    for i=1:nzcols
        idx = nzidx(i)-1;
        xx = mod(idx, length(Ws))+1;
        yy = mod(floor(idx/length(Ws)), length(Hs))+1;
        zz = floor(idx/length(Ws)/length(Hs))+1;
        fprintf(fp, '%g %g %g\n', Ws(xx), Hs(yy), Ds(zz)); % non-zero idx based on [W H D]
        
    end
    fclose(fp);
%     plot(squeeze(projected1(1,:,:)), squeeze(projected1(2,:,:)), '-');
    
end

