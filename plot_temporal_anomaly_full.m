% CASE = 4 % full annulus
% CASE = 4.1
% CASE = 5 % distribution
CASE = 6 % gmm results
% CASE_5_TYPE = -1  % anomaly
% CASE_5_TYPE = 0   % value
% CASE_5_TYPE = 1   % mahalanobis
% CASE_5_TYPE = 2   % mean
% CASE_5_TYPE = 3     % KLDiv
% CASE_5_TYPE = 4     % EMD

DISP_Y_RES = 11
SHOW_CASE = 'PASSAGE'
stepi = 1
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
    filepath ='/home/chenchu/volumes/TEMP/14.20_16.00_14.20x2_150523/tempdiff/'; TH=0.02
%     filepath = '/home/chenchu/volumes/TEMP/13.80_141219/tempdiff/' ; TH = 0.04
%     filepath = '/home/chenchu/volumes/TEMP/14.20_16.00/tempdiff/' ; TH = 0.02

    file_pattern = strcat('tempdiff_b%d_%d.raw')
    nfiles = 576
    period = 36
    W=151
    H=71
    D=56
    Ws = 11:2:81;
    Hs = 40:2:H;
    Ds = 1:2:D;
    starti=2
    endi = nfiles-period   
elseif CASE == 5         
%     filepath ='/home/chenchu/volumes/TEMP/14.20_16.00_14.20x2_150523/distribution_diff/';
%     filepath ='/home/chenchu/volumes/TEMP/14.20_16.00_14.20x2_150523/anomaly/';
%     filepath = '/home/chenchu/volumes/TEMP/13.80_141219/distribution_diff/' ;
    filepath = '/home/chenchu/volumes/TEMP/14.20_16.00/distribution_diff/' 
%     filepath = '/data/flow2/soumya/mean_13.8/anomaly/' 
    
    if CASE_5_TYPE==-1 % anomaly
        file_pattern = strcat('anomaly_b%d_%d.raw')
        TH = 2.8 %3.3296
        res = 1;
    elseif CASE_5_TYPE==0 % value
        file_pattern = strcat('distdiff_b%d_%d.raw')
        TH = 10
        res = 5;
    elseif CASE_5_TYPE==1 
        file_pattern = strcat('mahalanobis_b%d_%d.raw')
        TH = 5.5%3.3296
         res = 5;
    elseif CASE_5_TYPE==2
%         file_pattern = 'anomaly_b%d_%d.raw' %soumya
        file_pattern = 'mean_b%d_%d.raw'
        TH =  4 %3.3296
        res = 5;
    elseif CASE_5_TYPE==3
        file_pattern = 'kldiv_b%d_%d.raw' % KLDiv
        TH = 2;
        res = 5;
    elseif CASE_5_TYPE==4
        file_pattern = 'emd_b%d_%d.raw' % EMD
        TH = 3.33;
%         TH = 2;
        res = 5;
    end
    nfiles = 576
%     nfiles = 100
    period = 0
    W=151
    H=71
    D=56
    w1 = floor(W/res);
    h1 = floor(H/res);
    d1 = floor(D/res);
    Ws = 1:res:w1*res;
    Hs = 1:res:h1*res;
    Ds = 1:res:d1*res;
    starti=1
    endi = nfiles
elseif CASE == 6       
    %filepath =
    %'/home/chenchu/volumes/TEMP/uncertain/gmm_files_generated_pressure/anomaly/';     starti=6212    
%     filepath = '/home/chenchu/volumes/TEMP/uncertain/14.2x1/pressure_anomaly/out/';  starti = 21612
%     filepath = '/home/chenchu/volumes/TEMP/uncertain/14.2x2/pressure_anomaly/out/';  starti = 35012
%     filepath = '/home/chenchu/volumes/TEMP/uncertain/14.2x1/entropy_anomaly/raw/';  starti = 21612
    filepath = '/home/chenchu/volumes/TEMP/uncertain/14.2x2/entropy_anomaly/raw/';  starti = 35012
    file_pattern = 'emd_b%d_%d.raw'
    TH = 0.35 %3.3296
    res = 1;

    nfiles = 1339
    period = 360
    W=30
    H=14
    D=11
    D = ceil(D/DISP_Y_RES);
    Ws = 1:W;
    Hs = 1:H;
    Ds = 1:D;
    stepi = 10
    endi = starti+stepi*(nfiles-1)
end

DEGS = length(Ds)*36;
degmat = zeros(DEGS, length(starti:stepi:endi));
n = length(Ws)*length(Hs)*length(Ds);
for i=1:nfiles
    for b=1:36
        filename = sprintf(strcat(filepath, file_pattern), b, starti+stepi*(i-1))
        fp = fopen(filename, 'rb');
        data = fread(fp, n, 'float32');
        fclose(fp);

        list = find(abs(data)>TH);
        for j=1:length(list)
            idx = list(j)-1;
            xx = mod(idx, length(Ws))+1;
            yy = mod(floor(idx/length(Ws)), length(Hs))+1;            
            zz = floor(idx/length(Ws)/length(Hs))+1;
            zz = floor((zz-1)/DISP_Y_RES)+1; %!!!!!!!!!!
            y = zz+ length(Ds)*(b-1);
            
            if CASE==5 && (xx<4 || xx>24) %!!!!!!!!!!
                continue
            end
            if CASE==5 && CASE_5_TYPE==0
                degmat(y, i) = degmat(y, i) + data(list(j));    
            else
                degmat(y, i) = degmat(y, i) + 1;
%             degmat(y, i) = max(degmat(y, i), yy);
            end
        end
    end
    
end

figure
colormap(hot)
cmap = colormap(hot);
cmap = cmap([1:44,64],:);
cmap = flipud(cmap);
colormap(cmap)

% xrange = starti+period-1:endi+period-1; % temporal anomaly
xrange = 1:nfiles;
yrange = 1:length(Ds)*36;

maxval = max(max(degmat));
minval = min(min(degmat));
seg = (maxval-minval)/45;
h = imagesc(xrange, yrange, degmat, [minval-seg, (maxval)]);
yticklabel = 1:36;
xticklabel = starti:stepi:endi;
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
    'XColor', [.5 .5 .5],...
    'YColor', [.5 .5 .5] ,   ...
    'FontSize', 12, ...
    'GridLineStyle','-' ...
    )
xlabel('Time Step');
ylabel('Passage');


%% Draw 0-speed lines
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
out_image = sprintf('%s/anomaly.eps', filepath)
saveas(gca, out_image, 'psc2')
