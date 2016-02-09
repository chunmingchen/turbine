% data: (probe, time)
% time: (time)
function plot_probes(time, data, deglist)

DOWNSAMPLE=1

if DOWNSAMPLE
    time = time(1:25:end);
    data = data(:, 1:25:end);
    time = 1:length(time);
end
    
timesteps = size(data,2)
probes = size(data,1)
% [time, data]


% figure
max_stdev = max(std(data, 0, 1));

hold on
mag=40;
box on
for ii=1:probes
    %plot(time, data(ii,:)./3.0+ii, 'k')
    plot(time, data(ii,:)*mag+deglist(ii), 'k') %'Color', [.1, .1, .1])
end
% for t=1:timesteps
%     vals = data(:,t);
%     avg = mean(vals);
%     stdev = std(vals);
%     scatter(time(t)*ones(1,probes)', vals*mag + deglist' ,10*ones(1,probes)', abs(vals-avg)/stdev/1.96, 'filled');
% end
% c = colormap('hot');
% colormap( c(1:32, :))
hold off

if 1
    %ylim([0 360]);
    set(gca, 'YTickLabel',  [])
else
    offset = 0; %mean(data,2)*mag;
    set(gca, 'YTick', deglist+offset')
    if 1
        set(gca, 'YTickLabel',  deglist) %floor(linspace(0, 360, probes)))
    else
        set(gca, 'YTickLabel',  [])
    end
end
set(gca, 'FontSize', 20)
if DOWNSAMPLE
    xlabel('Time step (Down-sampled)')
else
    xlabel('Time step')
end
%set(gca, 'Color', [0 0 0])
ylabel('Normalized pressure offset by angle')
%saveas('pressure_probe_simulation.png')

    
end