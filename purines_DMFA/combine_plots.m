% This script reads the data for both GBM and normal cortex fluxes
% and combines the flux and confidence interval plots


%% Read the optimum fkux and the fux bounds

% 1 = GBM
% 2 = normal cortex

data_gbm = load('output_files\timestamp_20250303_1418\data.mat', ...
    'flux_cil', 'flux_ciu', 'flux_optim', 'rxns', 'tspan');
data_brain = load('output_files\timestamp_20250303_2209\data.mat', ...
    'flux_cil', 'flux_ciu', 'flux_optim');

flux_cil1 = data_gbm.flux_cil;
flux_ciu1 = data_gbm.flux_ciu;
flux_optim1 = data_gbm.flux_optim;

flux_cil2 = data_brain.flux_cil;
flux_ciu2 = data_brain.flux_ciu;
flux_optim2 = data_brain.flux_optim;

% common variables
rxns = data_gbm.rxns;
tspan = data_gbm.tspan;
clear data_gbm data_brain


%% Plot the results

folder = 'output_files\combine_plots';
pdf_folder = [folder, '\pdfs'];
mkdir(folder);
mkdir(pdf_folder);
for i = 1:numel(rxns)
    
    hold all

    patch([tspan fliplr(tspan)], [flux_cil2(i, :) fliplr(flux_ciu2(i, :))], ...
        validatecolor("#AFB0B0"), 'FaceAlpha', 0.35, 'EdgeColor', 'none')
    p1 = plot(tspan, flux_optim2(i, :), "Color", validatecolor("#5E686C"), ...
        "LineWidth", 1, 'DisplayName', 'Normal Cortex');
    patch([tspan fliplr(tspan)], [flux_cil2(i, :) fliplr(flux_ciu2(i, :))], ...
        validatecolor("#AFB0B0"), 'FaceAlpha', 0.35, 'EdgeColor', 'none')

    patch([tspan fliplr(tspan)], [flux_cil1(i, :) fliplr(flux_ciu1(i, :))], ...
        validatecolor("#FD866E"), 'FaceAlpha', 0.5, 'EdgeColor', 'none')
    p2 = plot(tspan, flux_optim1(i, :), "Color", validatecolor("#EE4928"), ...
        "LineWidth", 1, 'DisplayName', 'GBM');
    patch([tspan fliplr(tspan)], [flux_cil1(i, :) fliplr(flux_ciu1(i, :))], ...
        validatecolor("#FD866E"), 'FaceAlpha', 0.5, 'EdgeColor', 'none')

    hold off

    legend([p1 p2], 'FontSize', 15);
    xlim([0 tspan(end)]);
    ylim([0 inf]);
    title(rxns(i));
    xlabel('time (h)', 'FontSize', 18);
    ylabel('Flux (pmol/hr.mg-tissue)', 'FontSize', 18);
    
    filename = [folder, '\flux_CI', num2str(i), '.png'];
    saveas(gcf, filename);
    filename2 = [pdf_folder, '\flux', num2str(i), '.pdf'];
    saveas(gcf, filename2);
    close
      
end






