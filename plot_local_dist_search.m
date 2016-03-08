
CASE = 3.1

if CASE==3.1
    filepath = '/home/chenchu/volumes/TEMP/ld_search/TwoVariables_feature_11_neighbor_31/'
    file_pattern = 'search_result_14.2-13.8-turbine_%d.vti_Entropy_14.2-13.8-turbine_%d.vti_Pressure_31_0.9.raw'
    period = 144
    W=88
    H=509
    D=509
    Ws=1:W
    Hs=1:H
    Ds=1:H
    TH=0.16
    RADIAL_FRAME=true;
    istart = 300
    iend = 570
    istep = 10
    nfiles = (iend-istart)/istep+1;
end

if RADIAL_FRAME
    DEGS = 360;
else
    DEGS = length(Hs);
%     DEGS = length(Ds);
end
degmat = zeros(DEGS, nfiles);
n = length(Ws)*length(Hs)*length(Ds);
for i=1:nfiles
    filename = sprintf(strcat(filepath, file_pattern), istart+(i-1)*istep, istart+(i-1)*istep)
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
            degmat(yy+1, i) = degmat(yy+1, i) + 1;
%             degmat(zz+1, i) = degmat(zz+1, i) + 1;
        end
    end
    
end

figure
colormap(hot)
cmap = colormap(hot);
cmap = cmap([1:44,64],:);
cmap = flipud(cmap);
colormap(cmap)

maxval = max(max(degmat));
imagesc(degmat, [-maxval/45, maxval]);
xlabel('Time Step');
ylabel('Degree');
