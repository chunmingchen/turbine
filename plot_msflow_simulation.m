LARGE_FONT=0
filename = 'msflow.hs'
% filepath = '~/Dropbox/0turbine/141216_turbulent/';
% filepath = '/home/chenchu/Dropbox/0turbine/141219_simulation_result/oakley/'
% filepath = '/home/chenchu/Dropbox/0turbine/141219_simulation_result/glenn/'
% filepath = '/home/chenchu/Dropbox/0turbine/141222_simulation_result/' % 13.80
% filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_150112_turb_6201-20601/hs/'  %14.20
% filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_13.80_150127/'
% filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_14.00_150127_0-10585/'
% filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_150112_turb_continue_150216_20601-35001/hs/'
% filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_14.00_150127_continue_150212/'
% filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_150112_turb_continue_150216_20601-35001/hs/'
% filepath = '/home/chenchu/Dropbox/0turbine/150304/14.20/' % 14.20
% filepath = '/home/chenchu/Dropbox/0turbine/150304/16.00/'; %LARGE_FONT=1 % 16.00
% filepath = '/home/chenchu/Dropbox/0turbine/150427 simulation/14.00_16.00_13.80/' 
% filepath = '~/turbine_Stg/s35_noinj_14.20_16.00_13.80_150426/saved/hs/' %  16->13.80
filepath = '~/volumes/turbine_Stg_wd/s35_noinj_14.20_16.00_14.20_150426/saved/hs/' %  16->14.20
% filepath = '~/volumes/turbine_Stg_wd/s35_noinj_14.20_16.00_14.20x2_150523/saved/hs/'
% filepath = '~/turbine_Stg/s35_noinj_14.20_16.00_14.20x3_150528/saved/hs/' % 16->14.20x3
% filepath = '~/turbine_Stg/s35_noinj_14.20_16.00_14.20x4_150605/saved/hs/'
%filepath = '~/turbine_Stg/s35_noinj_14.20_16.00_14.20x5_150812/saved/hs/'% 16->14.20x5
%filepath = '~/Dropbox/0turbine/151209_single-passage/'
% filepath = '/home/chenchu/turbine_Stg_wd/single_passage/16.00_13.80/'
% filepath = '/home/chenchu/volumes/turbine_Stg_wd/single_passage/16.00_13.80x2/'
% filepath = '/home/chenchu/turbine_Stg_wd/single_passage/16.00_13.80x3/'
% filepath = '/home/chenchu/turbine_Stg_wd/single_passage/16.00_13.80x4/'
% filepath = '/home/chenchu/turbine_Stg_wd/single_passage/16.00_13.80x5/'
% filepath = '/home/chenchu/turbine_Stg_wd/single_passage/msflow_16.00_13.80_all.hs'; filename='';
% filepath = '/home/chenchu/volumes/TEMP/single_passage/16.00_13.80x3_12.00_10.00/'
% filepath = '/home/chenchu/volumes/TEMP/single_passage/16.00_13.80x3_12.00_11.00/'
% filepath = '/home/chenchu/volumes/TEMP/single_passage/16.00_13.80x3_12.00_10.50/'

rawdata = dlmread(strcat(filepath, filename));

% PAPER = 1
PAPER = 0
if LARGE_FONT
    fs = 30
else
    fs = 20
end

if PAPER
    rawdata = rawdata(1:25:end,:);
end

% first column: timestep
if PAPER 
    timestep = 1:size(rawdata,1);
elseif 1
    timestep = rawdata(:,1);
else
    timestep = (1:length(rawdata))/3600    
end
% 
data = rawdata(:,2:end); 

% plot(timestep, data, '--')
figure
if PAPER
    plot(timestep, data(:,3), 'k')  % 3: in flow 4: out flow
else
    plot(timestep, data)  % 3: in flow 4: out flow
end
%legend({'In flow', 'Out flow'})

% xlabel('revolution')
set(gca, 'FontSize', fs)
xlabel('Time step')
ylabel('Mass flow rate (kg/s)')
title('Mass Flow Rate')

if ~PAPER
    out_image = strcat(filepath, 'msflow.png');
    saveas(gca, out_image, 'png')
else
    out_image = strcat(filepath, 'msflow.eps');
    saveas(gca, out_image, 'psc2')
end
disp(sprintf('File saved to %s', out_image))
