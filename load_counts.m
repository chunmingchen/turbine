data=dlmread('count.txt');
% [time, block]
maxval = max(max(max(data)))
% 82 time steps, 36 blocks, 180 degrees
%figure('Color', [0.2 0.2 0.2])
% 
% figure
% colormap(hot);
% j=0
% for i=1:4:36
%     j=j+1;
%     subplot(3,3,j)    
%     %HeatMap(squeeze(data1(:,i,:))', 'LabelsWithMarkers', true, 'DisplayRange', 40, 'Colormap', redgreencmap(100, 'Interpolation', 'Sigmoid'))
%     plot(data(:,i)')
%     title(sprintf('block %d', i))
%     xlabel('Time Step')
%     ylabel('Count')
%     set(gca,'YDir','normal',  ...
%         'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
%         'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
%         'ytick', 0:10:100, ...
%         'color', 'w')
% %     ylim([15, 95])
% end
% saveas(gca, 'degrees.png', 'png');

j=0
for i=1:1:36
    j=j+1;
    subplot(6,6,j)    
    %HeatMap(squeeze(data1(:,i,:))', 'LabelsWithMarkers', true, 'DisplayRange', 40, 'Colormap', redgreencmap(100, 'Interpolation', 'Sigmoid'))
    plot(data(:,i)')
    title(sprintf('block %d', i))
    %xlabel('Time Step')
    ylabel('Count')
%     set(gca,'YDir','normal',  ...
%         'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
%         'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
%         'ytick', 70:10:100, ...
%         'color', 'w')
    ylim([0, maxval])
    %xlim([1, 576])
end
saveas(gca, 'count1.png', 'png');

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
%     ylim([0, maxval])
%     saveas(gca, sprintf('block %02d.png', i), 'png');
% %     pause
% end
% 
% [m,I] = max(data1, [], 3);
% I

