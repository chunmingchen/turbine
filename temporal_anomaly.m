filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_13.80_150127_regular_sampled/raw/'
filepattern = '14.2-13.8-turbine_%d.vti_Pressure_slice38.raw'
output_filepattern = 'tempdiff_slice38_%d.raw'
nfiles = 576
period = 144
W=509
H=509

%%
window = period;
first = true;
test = [];
test_fft = [];
test_dist = [];
data = zeros(W, H, window);
for i=1:(nfiles-window)
    if first
        for j=i:i+window-1
            filename = sprintf(strcat(filepath, filepattern), j)

            fp = fopen(filename, 'rb');
            data(:,:,j) = fread(fp, [W,H], 'float32');
            fclose(fp);
        end
    else
        data(:,:,1:end-1) = data(:,:,2:end);
        
        filename = sprintf(strcat(filepath, filepattern), i+window-1)
        fp = fopen(filename, 'rb');
        data(:,:,end) = fread(fp, [W,H], 'float32');
        fclose(fp);
    end
    test(i) = data(88,400,window);

    if 1 % debug
        mat = zeros(floor(window/2)+1, W,H);
        for x=1:W
            for y=1:H
%         x = 88 ; y = 400;
                s = squeeze(data(x,y,:));

                %display(sprintf('%d %d', ((t-1)*W+1) ,(t*W)))
                Y = fft(s);
                P2 = abs(Y)/window; % Fourier amplitude
                P1 = P2(1:window/2+1);
                P1(2:end-1) = 2*P1(2:end-1);

                mat(:,x,y) = P1;
            end
        end
        test_fft(:,i) = P1;
        if ~first
            % distance
            diff = mat-premat;
            dist = sqrt(sum(diff.^2,1));

            filename = sprintf(strcat(filepath, output_filepattern), i);
            fp = fopen(filename, 'wb');
            fwrite(fp, dist, 'float32');
            fclose(fp);
            
            test_dist(i) = dist(1,x,y);
        end
        premat = mat;
    end
    first = false;
end
