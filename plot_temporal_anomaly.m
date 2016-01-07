filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_13.80_150127_regular_sampled/raw/'
file_pattern = 'tempdiff_%d.raw'
nfiles = 576
period = 144
W=509
H=509
TH = 0.02

degmat = zeros(360, nfiles-period);
for i=2:nfiles-period
    filename = sprintf(file_pattern, i);
    fp = fopen(strcat(filepath, filename), 'rb');
    data = fread(fp, [W,H], 'float32');
    fclose(fp);
    
    [xlist,ylist] = find(data>TH);
    for j=1:length(xlist)
        x = xlist(j)-(W-1)/2;
        y = ylist(j)-(H-1)/2;
        deg = floor( (180/pi)*atan2(y,x) );
        if deg < 0
            deg = deg+360;
        end
        
        deg = deg + floor(i*360/144) - 50;
        deg = mod(deg, 360);
        degmat(deg+1, i) = degmat(deg+1, i)+1;
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
