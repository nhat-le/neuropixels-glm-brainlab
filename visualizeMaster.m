%% Visualize
files = dir('GLMResults\Forssmann_2017-11-01*ACA.mat');

%%
for i = 1:30
    fig = figure('Name', files(i).name);
    load(fullfile(files(i).folder, files(i).name))
    count = 1;
    lst = [1, 2, 3, 4, 5, 6, 7];
    fieldnames = fields(ws);
    for kCov = lst
        label = fieldnames{kCov};
        subplot(numel(lst), 1, count);
        %plot(ws.(label).tr, (ws.(label).data));
        errorbar(ws.(label).tr, ws.(label).data, sqrt(wvars.(label).data));
        hline(0)
        %ylim([-10 10])
        title(label);
        count = count + 1;
        ylim([-10 10])
    end
end