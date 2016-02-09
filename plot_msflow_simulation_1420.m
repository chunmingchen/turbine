filepath1 = '/data/flow2/turbine_Stg/s35_noinj_14.20_150112_turb_6201-20601/hs/'  %14.20
filepath2 = '/data/flow2/turbine_Stg/s35_noinj_14.20_150112_turb_continue_150216_20601-35001/hs/'
filepath3 = '/home/chenchu/Dropbox/0turbine/150304/14.20/' % 14.20
filepath = './' % output path
filename = 'msflow.hs'

rawdata1 = dlmread(strcat(filepath1, filename));
rawdata2 = dlmread(strcat(filepath2, filename));
rawdata3 = dlmread(strcat(filepath3, filename));
rawdata=[rawdata1;rawdata2;rawdata3];

PAPER = 1
fs = 20
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
title('Mass Flow Rate: In-flow')


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 4])
%set(gcf,'Papersize',[900 300])

if ~PAPER
    out_image = strcat(filepath, 'msflow.png');
    saveas(gca, out_image, 'png')
else
    out_image = strcat(filepath, 'msflow.eps');
    saveas(gca, out_image, 'psc2')
end
disp(sprintf('File saved to %s', out_image))
