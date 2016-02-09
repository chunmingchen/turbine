%target_path='/data/flow2/turbine_Stg/zDIR.P3D.rel.6201-11001/output_important/surfaces/'
target_path='./'
file=strcat(target_path, 'flowsepcount.txt');
data=dlmread(file);
data1 = reshape(data, 36, size(data,1)/36, size(data,2)); 
maxval = max(max(max(data1)))
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
%     imagesc(squeeze(data1(i,:,:))', [0 maxval])
%     title(sprintf('block %d', i))
%     xlabel('Time Step')
%     ylabel('Degree')
%     set(gca,'YDir','normal',  ...
%         'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
%         'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
%         'ytick', 0:10:100, ...
%         'color', 'w')
%     %ylim([15, 95])
% end
% %saveas(gca, 'flowsep.png', 'png');
% 
% j=0
% for i=1:1:36
%     j=j+1;
%     subplot(6,6,j)    
%     %HeatMap(squeeze(data1(:,i,:))', 'LabelsWithMarkers', true, 'DisplayRange', 40, 'Colormap', redgreencmap(100, 'Interpolation', 'Sigmoid'))
%     imagesc(squeeze(data1(i,:,:))', [0 maxval])
%     title(sprintf('Passage %d', i))
%     %xlabel('Time Step')
%     ylabel('Degree')
%     set(gca,'YDir','normal',  ...
%         'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
%         'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
%         'ytick', 70:10:100, ...
%         'color', 'w')
%     ylim([75, 95])
%     xlim([1, 30])
% end
% % saveas(gca, 'flowsep_large.png', 'png');

figure
j=0
timesteps = size(data1, 2)
degrees = size(data1,3)
for i=1:1:36
    j=j+1;
    subplot(6,6,j)    
    %HeatMap(squeeze(data1(:,i,:))', 'LabelsWithMarkers', true, 'DisplayRange', 40, 'Colormap', redgreencmap(100, 'Interpolation', 'Sigmoid'))
    %imagesc(squeeze(data1(i,:,:))', [0 maxval])
    hold on
    for t=1:timesteps
        for deg=1:degrees
            val = data1(i,t,deg);
            if val>0 %&& val>.1*maxval
                c=hsv2rgb([ max(0,min(1, deg/degrees)) ,1, 1 ]);  %
                scatter(t,deg,val/10,c,'fill')
            end
        end
    end
    hold off
    
    title(sprintf('Passage %d', i))
    %xlabel('Time Step')
    ylabel('Distance')
    set(gca,'YDir','normal',  ...
        'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
        'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
        ... %'ytick', 70:10:100, ...
        'color', 'w', ...
        'YDir', 'reverse')
%     set(gca, 'color', [0 0 0])
    ylim([1, degrees])
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

