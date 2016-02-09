% filepath = '~/Dropbox/0turbine/141216_turbulent/';
%filepath = '~/Downloads/';
% filepath = '/home/chenchu/Dropbox/0turbine/141219_simulation_result/oakley/'
% filepath = '/home/chenchu/Dropbox/0turbine/141219_simulation_result/glenn/'
% filepath = '/home/chenchu/Dropbox/0turbine/141222_simulation_result/' % 13.80
% filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_150112_turb_6201-20601/hs/' % 14.20
% filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_13.80_150127/'
% filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_14.00_150127/'
% filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_150112_turb_150212_continue/hs/'
% filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_150112_turb_continue_150216_20601-35001/hs/'
% filepath = '/home/chenchu/Dropbox/0turbine/150304/14.20/' % 14.20
% filepath = '/home/chenchu/Dropbox/0turbine/150304/16.00/' % 16.00
% filepath = '/media/chenchu/My Book/data/turbine_Stg/s35_noinj_14.20_16.00_13.80_150426/saved/hs/' %  16->13.80
% filepath = '/media/chenchu/My Book/data/turbine_Stg/s35_noinj_14.20_16.00_14.20_150426/saved/hs/' %  16->13.80
% filepath = '~/volumes/turbine_Stg_wd/s35_noinj_14.20_16.00_14.20_150426/saved/hs/' %  16->14.20
% filepath = '~/volumes/turbine_Stg_wd/s35_noinj_14.20_16.00_14.20x2_150523/saved/hs/' %  16->14.20x2
% filepath = '~/volumes/turbine_Stg_wd/s35_noinj_14.20_16.00_14.20x3_150528/saved/hs/' % 16->14.20x3
% filepath = '~/volumes/turbine_Stg_wd/s35_noinj_14.20_16.00_14.20x4_150605/saved/hs/'% 16->14.20x4
filepath = '~/volumes/turbine_Stg_wd/s35_noinj_14.20_16.00_14.20x5_150812/saved/hs/'% 16->14.20x5

name = 'probe_p_ring.1.hs'; CASE=2
%  name = 'probe_p_ring.2.hs'; CASE=3
%name = 'probe_p_ring.3.hs'; CASE=3
%name = 'probe_p_ring.4.hs'; CASE=1
%name = 'probe_p_ring.5.hs'; CASE=1

filename = strcat(filepath, name)

switch CASE
    case 1
        fin = fopen(filename, 'r');
        data0 = fscanf(fin, '%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ******** %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [36 Inf]);
        %                        1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0        1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  
        fclose(fin);
        data0 = data0';
        size(data0)

        time = data0(:,1)';
        probes = size(data0, 2)-1
        data = data0(:,2:probes+1)'; % data(probe, time)
        for i=1:probes3600
            data(i,:)=data(i,:)-(i-1);
        end

        deg = 0:10:359;
        if 0
            deg(21)=[]
        end

        figure
        plot_probes(time, data, deg)

        plot_fft(time,data(1,:));

    case 2
        data0 = dlmread(filename);
        size(data0)

        time = data0(:,1)';
        probes = size(data0, 2)-1
        data = data0(:,2:probes+1)'; % data(probe, time)
        for i=1:probes
            data(i,:)=data(i,:)-(i-1);
        end

        figure
        plot_probes(time, data, [10, 70, 100, 160, 190, 250, 280, 340]);
        plot_fft(time, data(2,:));

    case 3 % probe 3
        data0 = dlmread(filename);

        size(data0)

        time = data0(:,1)';
        probes = size(data0, 2)-1
        data = data0(:,2:probes+1)'; % data(probe, time)
        for i=1:probes
            data(i,:)=data(i,:)-(i-1);
        end

        figure
        plot_probes(time, data, [0 90 180 270]);
        plot_fft(time, data(2,:));

    otherwise
        disp 'CASE unknown'
end 
saveas(gca, strcat(filename, '.png'))
saveas(gca, strcat(filename, '.eps'), 'psc2')
    
    % 
%     timesteps = size(data0,1)
%     probes = size(data0,2)-1
%     
%     [time, data]
%     for i=1:probes
%         data(:,i)=data(:,i)-(i-1);
%     end
%     
%     figure
%     max_stdev = max(std(data, 0, 2));
%     deg = 360/probes;
%     
%     hold on
%     for t=1:timesteps
%         vals = data(t, :);
%         avg = mean(vals);
%         stdev = std(vals);
%         scatter(time(t)*ones(1,probes), vals/3 + (1:probes)-1 ,10*ones(1,probes), abs(vals-avg)/stdev, 'filled');
%     end
%     hold off
%     c = colormap('hot');
%     colormap( c(1:20, :))
%     xlabel('Time Step')
%     set(gca, 'YTick', [0, 1, 2, 3, 4, 5, 6, 7])
%     set(gca, 'YTickLabel', [0, 45, 90, 135, 180, 225, 270, 315])
%     ylabel('Degree')
%     saveas('pressure_probe_simulation.png')
%     
%     
      