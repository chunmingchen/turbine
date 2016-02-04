CASE = 4 % full annulus
CASE = 4.1
SHOW_CASE = 'PASSAGE'
if CASE==4
    filepath = '/media/chenchu/TEMP/14.20_13.80_single/tempdiff/'
    file_pattern = strcat('tempdiff_b%d_%d.raw')
    nfiles = 534
    period = 144
    W=151
    H=71
    D=56
    Ws = 11:2:81;
    Hs = 40:2:H;
    Ds = 1:2:D;
    TH = 0.04
    starti=2
    endi = nfiles-period
elseif CASE == 4.1        
    filepath = '/home/chenchu/volumes/TEMP/14.20_16.00_14.20x2_150523/tempdiff/'
    file_pattern = strcat('tempdiff_b%d_%d.raw')
    nfiles = 576
    period = 36
    W=151
    H=71
    D=56
    Ws = 11:2:81;
    Hs = 40:2:H;
    Ds = 1:2:D;
    TH = 0.02
    starti=2
    endi = nfiles-period        
end
    
DEGS = length(Ws)*36;
degmat = zeros(DEGS, endi)-1;
n = length(Ws)*length(Hs)*length(Ds);
for i=starti:endi
    for b=1:36
        filename = sprintf(strcat(filepath, file_pattern), b, i)
        fp = fopen(filename, 'rb');
        data = fread(fp, n, 'float32');
        fclose(fp);

        list = find(data>TH);
        for j=1:length(list)
            idx = list(j)-1;
            xx = mod(idx, length(Ws));
            yy = mod(floor(idx/length(Ws)), length(Hs));
            zz = floor(idx/length(Ws)/length(Hs));
            y = xx+1+ length(Ws)*(b-1);
            degmat(y, i) = degmat(y, i) + 1;
%             degmat(y, i) = max(degmat(y, i), yy);
        end
    end
    
end

figure
colormap(hot)
cmap = colormap(hot);
cmap = cmap([1:44,64],:);
cmap = flipud(cmap);
colormap(cmap)

xrange = starti+period-1:endi+period-1;
yrange = 1:length(Ds)*36;
imagesc(xrange, yrange, degmat);
yticklabel = 1:36;
% set(gca,'YTickLabel',{'9' '8' '7' '6' '5' '4' '3' '2' '1' '0'})
set(gca, ...
            'ytick', 1:length(Ds):(36*length(Ds)-1), ...
            'yticklabel', yticklabel, ...
            'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
            'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
            'GridLineStyle','-' ...
        ) 
xlabel('Time Step');
ylabel('Passage');


hold on

for i=floor(starti/144):ceil(endi/144)
    if strcmp(SHOW_CASE, 'PASSAGE')
        x=[(i-1)*144+1; i*144+1];
        y=[1; 36*length(Ds)+1];
        plot([x(1);x(2)], [y(1);y(2)], ':', 'color', [.5 .5 .5])
    end
end

hold off

cbr = colorbar;
% set(cbr, 'YTick', [0])