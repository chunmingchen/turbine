p=0

fname = sprintf('ring%d.raw', p)
fp = fopen(fname, 'r')
a=fread(fp, inf, 'float32');
fclose(fp);
data =reshape(a, 56*36, length(a)/56/36);

minval = min(a);

ts = size(data, 2);
% for t=1:ts
%     for k=0:55
%         group = data( (0:35)*56+k+1, t );
% 
%         if 1
%             % anomaly analysis
% %             m=mean(group);
% %             s=std(group);
% %             ai=find(abs(group-m)/s < 2.3);
% %             if length(ai)>0
% %                 %data((ai-1)*56+k+1,t)= minval;
% %                 %data((ai-1)*56+k+1,t)= 0;%data((ai-1)*56+k+1,t)+1;
%             end
%         else
%             m=mean(group);
%             s=std(group);
%             data( (0:35)*56+k+1, t ) = (group-m)/s;
%         end
%     end
% end
% for k=1:size(data,1)
%     group = data(k, :);
%     m = mean(group);
%     s = std(group);
%     data(k, :) = ((group-m)/s);
% end


figure
if 0
    cmap = colormap(hot);
    cmap = cmap([1:44,64],:);
    cmap = flipud(cmap);
    colormap(cmap)
elseif 0
    colormap (hot)
%         c = colormap 
%         c= flipud(c)
%         colormap (c)
else
    %colormap (winter)
    colormap (parula)
end
imagesc(data)
colorbar
yticklabel = 1:36;
ylabel('Passage')
set(gca, 'ytick', 1:56:(36*56-1), ...
        'yticklabel', yticklabel)
title(fname)

hold on
for i=1:4
    x=[(i-1)*144+1; i*144+1];
    y=[1; 36*56+1];
    plot([x(1);x(2)], [y(1);y(2)], '--', 'color', [.5 .5 .5])
end
hold off

% plot_probe(