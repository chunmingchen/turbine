i=1
close all
dataall=zeros(36, 50);
for b=1:1:36
    %subplot(9,1,i)
    %i=i+1
    figure
    filename=sprintf('probed_pressure/PressureProbeBlock%d.txt', b)
    data=dlmread(filename);
    % [time, block]
    maxval = max(max(max(data)))
    pos = data(2:end,1:3);
    val = data(2:end,4:end);
    x=1:50
    id = 1:8:64
    plot(x, val(id, x)'+ones(50,1)*[1:8])
    %ylim([0, 2])
    title(filename)
    %saveas(gca, sprintf('probed_pressure/pressure1_%02d.png', b))
    id = 1
    dataall(b, :) = val(id, :);
end

i=1
close all
figure('Color',[0 0 0]);
for id=20:20
    data_all=zeros(36, 192);
    pos_all = [];
    for b=1:36
        %i=i+1
        filename=sprintf('probed_pressure_New1/PressureProbeBlock%d.txt', b)
        data=dlmread(filename);
        % [time, block]
        %maxval = max(max(max(data)))
        pos = data(2:end,1:3);
        val = data(2:end,4:end);
        
        data_all(b, :) = val(id, :);
        pos_all = pos(id, :);
    end
    if 0
        x=1:192;
        plot(x, data_all(1:1:36,:)'+ones(192,1)*[1:36])
        title(sprintf('all time steps at pos %f %f %f', pos_all(1), pos_all(2), pos_all(3)))
        xlabel('Time step')
        ylabel('Block')
    %     saveas(gca, sprintf('probed_pressure_New1/pressure_id_%02d.png', id))
    else
        %!!!!!!! new 
        data_all1 = data_all(1:36,:);
        probes = size(data_all1,1);
        timesteps = size(data_all1,2);
        data_all2 = zeros(probes, timesteps*25);
        for i=1:probes
            data_all2(i,:) = interp1(0:25:(timesteps-1)*25, data_all1(i,:), 0:(timesteps*25-1));
        end
        plot_probes(6201+(0:2:timesteps*25-1), data_all2(:,1:2:timesteps*25))
    end
end

if 0

    b=28
    %i=i+1
    figure
    filename=sprintf('probed_pressure/PressureProbeBlockStable%d.txt', b)
    filename=sprintf('probed_pressure/PressureProbeBlock%dall.txt', b)
    data=dlmread(filename);
    % [time, block]
    maxval = max(max(max(data)))
    pos = data(2:end,1:3);
    val = data(2:end,4:end);
    x=1:192
    id = 1:4:64;
    plot(x, val(id, x)'+ones(192,1)*[1:16])
    plot(x, val(:, x)')
    %ylim([0, 2])
    title(filename)
    %saveas(gca, sprintf('probed_pressure/pressure1_%02d.png', b))
    id = 5
    dataall(b, :) = val(id, :);
end
% 
% % 82 time steps, 36 blocks, 180 degrees
% %figure('Color', [0.2 0.2 0.2])
% % 
% % figure
% % colormap(hot);
% % j=0
% % for i=1:4:36
% %     j=j+1;
% %     subplot(3,3,j)    
% %     %HeatMap(squeeze(data1(:,i,:))', 'LabelsWithMarkers', true, 'DisplayRange', 40, 'Colormap', redgreencmap(100, 'Interpolation', 'Sigmoid'))
% %     plot(data(:,i)')
% %     title(sprintf('block %d', i))
% %     xlabel('Time Step')
% %     ylabel('Count')
% %     set(gca,'YDir','normal',  ...
% %         'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
% %         'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
% %         'ytick', 0:10:100, ...
% %         'color', 'w')
% % %     ylim([15, 95])
% % end
% % saveas(gca, 'degrees.png', 'png');
% 
% j=0
% for i=1:1:36
%     j=j+1;
%     subplot(6,6,j)    
%     %HeatMap(squeeze(data1(:,i,:))', 'LabelsWithMarkers', true, 'DisplayRange', 40, 'Colormap', redgreencmap(100, 'Interpolation', 'Sigmoid'))
%     plot(data(:,i)')
%     title(sprintf('block %d', i))
%     %xlabel('Time Step')
%     ylabel('Count')
% %     set(gca,'YDir','normal',  ...
% %         'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
% %         'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
% %         'ytick', 70:10:100, ...
% %         'color', 'w')
%     ylim([0, maxval])
%     xlim([1, 20])
% end
% saveas(gca, 'count1.png', 'png');
% 
% % 
% % %colorbar
% % 
% % figure
% % colormap(hot);
% % for i=1:36
% %     %HeatMap(squeeze(data1(:,i,:))', 'LabelsWithMarkers', true, 'DisplayRange', 40, 'Colormap', redgreencmap(100, 'Interpolation', 'Sigmoid'))
% %     imagesc(squeeze(data1(i,:,:))', [0 maxval])
% %     title(sprintf('block %d', i))
% %     xlabel('Time Step')
% %     ylabel('Degree')
% %     set(gca,'YDir','normal',  ...
% %         'xgrid', 'on', 'xcolor',[.2 .2 .2], ...
% %         'ygrid', 'on', 'ycolor',[.2 .2 .2], ...
% %         'ytick', 0:10:180, ...
% %         'color', 'w')
% %     ylim([0, maxval])
% %     saveas(gca, sprintf('block %02d.png', i), 'png');
% % %     pause
% % end
% % 
% % [m,I] = max(data1, [], 3);
% % I
% 
