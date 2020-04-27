% To visualize results of GLM fitting for multiple cells

% Load data
load ws_wvar_VISp_Cori1214ACA.mat
vars = fieldnames(ws_all);
nCovar = numel(vars);

lst = [1, 2, 3, 4, 5, 6, 7];

for i = 1:numel(ws_all)
    fig = figure(1); clf;
    count = 1;
    for kCov = lst
        label = vars{kCov};
        subplot(numel(lst), 1, count);
        %plot(ws.(label).tr, (ws.(label).data));
        errorbar(ws_all(i).(label).tr, ws_all(i).(label).data, sqrt(wvars_all(i).(label).data));
        hline(0)
        %ylim([-10 10])
        title([label num2str(i)]);
        count = count + 1;
        ylim([-10 10])
    end
    figname = sprintf('ACApcell%d.fig', i);
    savefig(gcf, figname, 'compact')
end