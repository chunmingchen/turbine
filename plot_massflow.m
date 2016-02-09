% dir = '/data/flow2/turbine_Stg/s35_noinj_13.80_141219_turb_6201-20601/output_important/'
% dir = '/data/flow2/turbine_Stg/s35_noinj_14.20_13.80_150127/'


% filepath = '/home/chenchu/Dropbox/0turbine/141222_simulation_result/' % 13.80
% filepath = '/data/flow2/turbine_Stg/s35_noinj_13.80_150112_turb_nostall/hs/'
%filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_13.80_150127/'
%filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_14.00_150127/'
% filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_150112_turb_150212_continue/hs/'
% filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_14.00_150127_continue_150212/'
% filepath = '/data/flow2/turbine_Stg/s35_noinj_14.20_150112_turb_continue_150216_20601-35001/hs/'
% filepath = '/home/chenchu/Dropbox/0turbine/150304/14.20/' % 14.20
%filepath = '/home/chenchu/Dropbox/0turbine/150304/16.00/' % 16.00
%dir = '/media/chenchu/My Book/data/turbine_Stg/s35_noinj_14.20_16.00_13.80_150426/saved/' %  16->13.80
% dir = '~/turbine_Stg/s35_noinj_14.20_16.00_14.20_150426/saved/' %  16->14.20
dir = '~/turbine_Stg/s35_noinj_14.20_16.00_14.20x3_150528/saved/' % 16->14.20x3

if 0 % paper:
    x = load(strcat(dir, 'inflow.txt'));
    y = load(strcat(dir, 'outflow.txt'));
    % x = x*106.22;
    % y = y*106.22;
    fs = 20

    figure
    plot(x,'k');
    ti=title('Mass Flow Rate: In-flow');
    xl=xlabel('Time step (Down-sampled)');
    yl=ylabel('Mass flow rate (kg/s)');
    set(gca, 'FontSize', fs)
    set(xl, 'FontSize', fs)
    set(yl, 'FontSize', fs)
    set(ti, 'FontSize', fs)

    out_image = strcat(dir , 'mass_flow.eps'); %strcat(filepath, 'msflow.png');
    disp(sprintf('File saved to %s', out_image))
    saveas(gca, out_image, 'psc2')

else
    %%
    x = load(strcat(dir, 'inflow.txt'));
    y = load(strcat(dir, 'outflow.txt'));
    % x = x*106.22;
    % y = y*106.22;

    plot(x,'r');
    hold on
    plot(y,'b');
    hold on
    legend('in-flow','out-flow');
    title('Mass Flow Rate: In-flow and Out-flow');
    xlabel('time steps');
    ylabel('mass flow rate (kg/s)');


end