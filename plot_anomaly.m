timesteps = 576 % 193 %size(data1, 2)
% target_path='/data/flow2/turbine_Stg/zDIR.P3D.rel.6201-11001/output/'
target_path='/data/flow2/turbine_Stg/s35_noinj_13.80_141219_turb_6201-20601/output_important/anomaly/'
var = 'Pressure'
% var = 'Entropy'
% var = 'TotalPressure'
% var = 'Density'
% var = 'Temperature'
% var = 'Velocity'
data1=zeros(36*56, 16, timesteps);
for t=1:timesteps
    file=sprintf('%s/Salient_%s_iso3.5_%d.txt', target_path, var, t-1);
    data=dlmread(file);
    data1(:,:,t)=data;
end

% remove small area
conn = 6
P = 500
P1 = 20000
% P1= 3406207*.5  %entropy
% P1 = 19885237 *.95  %velocity
% P1 = 273366 *.05 %density
if 0
    BW2 = bwareaopen(data1, P, conn);
    data1 = data1.*BW2;
elseif 1
    % Determine the connected components:
    CC = bwconncomp(data1, conn);
    % Compute the area of each component:
    if 0
        S = regionprops(CC, 'Area');
        % Remove small objects:
        L = labelmatrix(CC);
        BW2 = ismember(L, find([S.Area] >= P));
        data1 = data1.*BW2
    else
        S = regionprops(CC, 'Area');
        S1 = regionprops(CC, data1, 'MeanIntensity');
        % Remove small objects:
        L = labelmatrix(CC);
        maxsum = max([S.Area].*[S1.MeanIntensity])
        BW2 = ismember(L, find(([S.Area].*[S1.MeanIntensity]) >= P1));
        data1 = data1.*BW2;
    end
end


figure
degrees = size(data1,1);
% colormap(hot)

        cmap = colormap(hot);
        cmap = cmap([1:44,64],:);
        cmap = flipud(cmap);
        colormap(cmap)
        
data2 = squeeze(sum(data1,2));
maxval = max(max(max(data2)))
imagesc(sqrt(data2), [0 sqrt(maxval)])
% imagesc((data2), [0 (maxval)])
%xlabel('Time Step')
ylabel('Passage')
set(gca,'YDir','normal',  ...
    'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
    'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
    'ytick', 1:56:(36*56-1), ...
    'yticklabel', 1:36, ...
    'xtick', 1:5:timesteps, ...
    'color', 'w', ...
    'YDir', 'reverse', ...
    'XColor', [.6 .6 .6],...
    'YColor', [.6 .6 .6]    ...
    )
set(gca, 'color', [0 0 0])
ylim([1, degrees])
xlim([1, timesteps])

if 0
    maxval = max(max(max(data1)))
    % 193 time steps, 36 x 56 degrees, 10 segs

    % merge
    if 1
        out_segs = 8
        segs = size(data1,2);
        data2=zeros(size(data1,1), out_segs, size(data1,3));
        step = ceil(segs/out_segs);
        for i=1:out_segs
            merge = ((i-1)*step+1):min(segs,(i*step))
            data2(:,i,:) = sum(data1(:,merge,:),2);
        end 
        data1 = data2;
    end

    figure
    j=0;
    segs = size(data1,2)
    degrees = size(data1,1)
    if 0
        for i=1:1:segs
            j=j+1;
            %HeatMap(squeeze(data1(:,i,:))', 'LabelsWithMarkers', true, 'DisplayRange', 40, 'Colormap', redgreencmap(100, 'Interpolation', 'Sigmoid'))
            %imagesc(squeeze(data1(i,:,:))', [0 maxval])
            data = squeeze(data1(1:degrees, i, 1:timesteps));
            [degs,ts,val] = find(data);
            scatter(ts,degs,(val/maxval).^2*200,  (i/segs)*ones(length(ts),1),'fill')

            title(sprintf('X %d', i))
            %xlabel('Time Step')
            ylabel('Passage')
            set(gca,'YDir','normal',  ...
                'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
                'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
                'ytick', 1:56:(36*56-1), ...
                'yticklabel', 1:36, ...
                'color', 'w', ...
                'YDir', 'reverse')
            set(gca, 'color', [0 0 0])
            ylim([1, degrees])
            xlim([1, timesteps])
            hold on
        end
        hold off
    else
    %     cmap = colormap(hot);
    %     cmap = cmap(8:56,:);
    %     cmap = flipud(cmap);
        colormap(hot);
        for i=1:1:segs
            j=j+1;
            subplot(2,ceil(segs/2),j)    
            %HeatMap(squeeze(data1(:,i,:))', 'LabelsWithMarkers', true, 'DisplayRange', 40, 'Colormap', redgreencmap(100, 'Interpolation', 'Sigmoid'))
            %imagesc(squeeze(data1(i,:,:))', [0 maxval])
            data = squeeze(data1(1:degrees, i, 1:timesteps));
            if 0
                [degs,ts,val] = find(data);
                scatter(ts,degs,(val/maxval)*200,  sqrt(val/maxval),'fill')
            else
                imagesc(sqrt(data), [0 sqrt(maxval)])
            end

            title(sprintf('%.2f%%-%.2f%% chord from leading edge', (i-ceil(segs/2)-1)/3, (i-ceil(segs/2))/3))
            %xlabel('Time Step')
            ylabel('Passage')
            set(gca,'YDir','normal',  ...
                'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
                'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
                'ytick', 1:(56*4):(36*56-1), ...
                'yticklabel', 1:4:36, ...
                'color', 'b', ...
                'YDir', 'reverse')
            ... %set(gca, 'color', [0 0 0])
            ylim([1, degrees])
            xlim([1, timesteps])
        end
    end
    % 
    % %colorbar
    % 
    % figure
    % colormap(hot);
    % for i=1:36
    %     %HeatMap(squeeze(data1(:,i,:))', 'LabelsWithMarkers', true, 'DisplayRange', 40, 'Colormap', redgreencmap(100, 'Interpolation', 'Sigmoid'))
    %     imagesc(squeeze(data1(i,:,:))', [0 maxval])
    %     title(sprintf('block %d', i))
    %     xlabel('Time Step')
    %     ylabel('Degree')
    %     set(gca,'YDir','normal',  ...
    %         'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
    %         'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
    %         'ytick', 0:10:180, ...
    %         'color', 'w')
    %     ylim([0, 95])
    %     saveas(gca, sprintf('block %02d.png', i), 'png');
    % %     pause
    % end
    % 
    % [m,I] = max(data1, [], 3);
    % I

end