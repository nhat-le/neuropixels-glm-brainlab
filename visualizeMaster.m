%% Visualize
files = dir('GLMResults\Radnitz_2017-01-08_unit626_ACA.mat');

%%
for i = 1:50
    fig = figure('Name', files(i).name);
    load(fullfile(files(i).folder, files(i).name))
    count = 1;
    lst = [1, 3, 4, 5, 7];
    fieldnames = fields(ws);
    for kCov = lst
        label = fieldnames{kCov};
        subplot(numel(lst), 1, count);
        %plot(ws.(label).tr, (ws.(label).data));
        errorbar(ws.(label).tr(1:10:end), ws.(label).data(1:10:end), sqrt(wvars.(label).data(1:10:end)));
        hline(0)
        %ylim([-10 10])
        title(label);
        count = count + 1;
        ylim([-5 5])
    end
end