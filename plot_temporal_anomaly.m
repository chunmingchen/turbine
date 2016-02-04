% CASE=1
% CASE=2 % single passage
% CASE = 3 % anomaly full
% CASE = 1.5 % single passage all points
% CASE = 4 % single block from full
CASE = 4.1 % full annulus

if CASE==1
    filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_13.80_150127_regular_sampled/raw/'
    file_pattern = 'tempdiff_%d.raw'
    nfiles = 576
    period = 144
    W=1
    H=509
    D=509 
    TH = 0.01
    RADIAL_FRAME=true
    starti=2
    endi = nfiles-period
elseif CASE==1.5
    %filepath = '/media/chenchu/TEMP/14.20_13.80_raw_low/tempdiff/'
    %file_pattern = 'tempdiff_%d.raw'
    %TH = 0.02
%     starti=2
    filepath = '/media/chenchu/TEMP/14.20_13.80_raw_low/li2016/'
    file_pattern = 'li2016_%d.raw'
    TH = -0.96
    starti=1
    
    nfiles = 534
    period = 144
    W=35
    H=204
    D=204
    RADIAL_FRAME=true
    endi = nfiles-period*2
elseif CASE==2
    slice='20' % '40'
    filepath = '//media/chenchu/TEMP/14.20_13.80_single/'
%     file_pattern = strcat('tempdiff_slice',slice,'_%d.raw')
    file_pattern = strcat('tempdiff_b10_slice',slice,'_%d.raw')
    nfiles = 534
    period = 144
    W=1
    H=71
    D=56
    TH = 0.02
    RADIAL_FRAME=false
    starti=2
    endi = nfiles-period
elseif CASE==3
    filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_13.80_150127_regular_sampled/raw/'
    file_pattern = 'anomaly/slice38_%d.raw'
    nfiles = 575
    period = 144
    W=1
    H=509
    D=509
    TH = 2.5%3.3296
    RADIAL_FRAME=true
    starti=1
    endi = nfiles
elseif CASE==4 
    block = '10';
    filepath = '/media/chenchu/TEMP/14.20_13.80_single/tempdiff/'
    file_pattern = strcat('tempdiff_b', block , '_%d.raw')
    nfiles = 534
    period = 144
    W=151
    H=71
    D=56
    Ws = 11:2:81;
    Hs = 40:2:H;
    Ds = 1:2:D;
    TH = 0.08
    RADIAL_FRAME=false
    starti=2
    endi = nfiles-period
elseif CASE==4.1    
    block = '10';
    filepath = '/home/chenchu/volumes/TEMP/14.20_16.00_14.20x2_150523/tempdiff/'
    file_pattern = strcat('tempdiff_b', block , '_%d.raw')
    nfiles = 576
    period = 144
    W=151
    H=71
    D=56
    Ws = 11:2:81;
    Hs = 40:2:H;
    Ds = 1:2:D;
    TH = 0.02
    starti=2
    endi = nfiles-period    
    RADIAL_FRAME=false    
end
    
if RADIAL_FRAME
    DEGS = 360;
else
    DEGS = length(Ws);
end
degmat = zeros(DEGS, endi);
n = length(Ws)*length(Hs)*length(Ds);
for i=starti:endi
    filename = sprintf(strcat(filepath, file_pattern), i)
    fp = fopen(filename, 'rb');
    data = fread(fp, n, 'float32');
    fclose(fp);
    
    list = find(data>TH);
    for j=1:length(list)
        idx = list(j)-1;
        xx = mod(idx, length(Ws));
        yy = mod(floor(idx/length(Ws)), length(Hs));
        zz = floor(idx/length(Ws)/length(Hs));
        if RADIAL_FRAME
            x = yy-(length(Hs)-1)/2;
            y = zz-(length(Ds)-1)/2;
            deg = floor( (180/pi)*atan2(y,x) );
            if deg < 0
                deg = deg+360;
            end

            deg = deg + floor(i*360/period) - 50;
            deg = mod(deg, 360);
            degmat(deg+1, i) = degmat(deg+1, i)+1;
        else
            degmat(xx+1, i) = degmat(xx+1, i) + 1;
        end
    end
    
end

figure
colormap(hot)
cmap = colormap(hot);
cmap = cmap([1:44,64],:);
cmap = flipud(cmap);
colormap(cmap)

imagesc(degmat);
xlabel('Time Step');
ylabel('Degree');
