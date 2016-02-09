data=dlmread('degrees.txt');
data1 = reshape(data, 36, size(data,1)/36, 180); 
maxval = max(max(max(data1)))
% 82 time steps, 36 blocks, 180 degrees
%figure('Color', [0.2 0.2 0.2])

figure
colormap(hot);
j=0
for i=1:4:36
    j=j+1;
    subplot(3,3,j)    
    %HeatMap(squeeze(data1(:,i,:))', 'LabelsWithMarkers', true, 'DisplayRange', 40, 'Colormap', redgreencmap(100, 'Interpolation', 'Sigmoid'))
    imagesc(squeeze(data1(i,:,:))', [0 maxval])
    title(sprintf('block %d', i))
    xlabel('Time Step')
    ylabel('Degree')
    set(gca,'YDir','normal',  ...
        'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
        'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
        'ytick', 0:10:100, ...
        'color', 'w')
    ylim([15, 95])
end
saveas(gca, 'degrees.png', 'png');

j=0
for i=1:1:36
    j=j+1;
    subplot(6,6,j)    
    %HeatMap(squeeze(data1(:,i,:))', 'LabelsWithMarkers', true, 'DisplayRange', 40, 'Colormap', redgreencmap(100, 'Interpolation', 'Sigmoid'))
    imagesc(squeeze(data1(i,:,:))', [0 maxval])
    title(sprintf('Passage %d', i))
    %xlabel('Time Step')
    ylabel('Degree')
    set(gca,'YDir','normal',  ...
        'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
        'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
        'ytick', 70:10:100, ...
        'color', 'w')
    ylim([75, 95])
    xlim([1, 500])
end
saveas(gca, 'degrees_large.png', 'png');

figure
j=0
timesteps = size(data1, 2)
degrees = size(data1,3)
for i=1:1:36
    j=j+1;
    subplot(6,6,j)    
    %HeatMap(squeeze(data1(:,i,:))', 'LabelsWithMarkers', true, 'DisplayRange', 40, 'Colormap', redgreencmap(100, 'Interpolation', 'Sigmoid'))
    %imagesc(squeeze(data1(i,:,:))', [0 maxval])
    if 0
        hold on
        for t=1:timesteps
            for deg=1:degrees
                val = data1(i,t,deg);
                if val>0 %&& val>.1*maxval
                    c=hsv2rgb([ max(0,min(1, 1-(deg-70)/15)) ,1, 1 ]);  %
                    %scatter(t,deg,(val*30000),c,'fill')
                    scatter(t,deg,(val/maxval*1500),c,'fill')
                end
            end
        end
        hold off
    else
        [t,deg,vals] = find(squeeze(data1(i,:,:)));
        n= length(t);
        c = zeros(n,3);
        for x=1:n
            c(x,:)=hsv2rgb([ max(0,min(1, 1-(deg(x)-70)/15)) ,1, 1 ]);  %
        end
        scatter(t,deg,250*(vals/maxval),c,'fill')
    end
    
    title(sprintf('Passage %d', i))
    %xlabel('Time Step')
    ylabel('Degree')
    set(gca,'YDir','normal',  ...
        'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
        'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
        'ytick', 70:10:100, ...
        'color', 'w')
    ylim([75, 95])
%     set(gca, 'color', [0 0 0])
    xlim([1, timesteps])
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

