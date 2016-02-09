filepath = {'~/turbine_Stg_wd/s35_noinj_14.20_16.00_14.20_150426/saved/hs/', ... %  16->14.20
            '~/turbine_Stg_wd/s35_noinj_14.20_16.00_14.20x2_150523/saved/hs/', ... %  16->14.20x2
            '~/turbine_Stg_wd/s35_noinj_14.20_16.00_14.20x3_150528/saved/hs/', ... % 16->14.20x3
            '~/turbine_Stg_wd/s35_noinj_14.20_16.00_14.20x4_150605/saved/hs/', ...% 16->14.20x4
            '~/turbine_Stg_wd/s35_noinj_14.20_16.00_14.20x5_150812/saved/hs/'}% 16->14.20x5

name = 'probe_p_ring.1.hs'; CASE=2
%  name = 'probe_p_ring.2.hs'; CASE=3

time_list = cell(5,1);
data_list = cell(5,1);
probe = 2;

for j = 1:5
    filename = strcat(filepath{j}, name)

    data0 = dlmread(filename);

    size(data0)

    time = data0(:,1)';
    probes = size(data0, 2)-1
    data = data0(:,2:probes+1)'; % data(probe, time)
    for i=1:probes
        data(i,:)=data(i,:)-(i-1);
    end

    time_list{j} = time;
    data_list{j} = data(probe,:);

end

time_list
time = [time_list{1} time_list{2} time_list{3} time_list{4} time_list{5}];
data = [data_list{1} data_list{2} data_list{3} data_list{4} data_list{5}];

plot(time/3600, data);
% plot_fft(time, data);
plot_windowed_fft
% saveas(gca, sprintf('wfft_probe2_%d.png', probe) )
% plot_morlet