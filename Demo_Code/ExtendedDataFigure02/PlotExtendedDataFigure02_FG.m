clear; close all;

%% update your main directory to manuscript folder 
maindir = 'Y:\labs\singer\Xiao\Code\projects\nuri_manuscript_figs\ManuscriptCode\InhibitoryGating2023';
addpath(genpath(fullfile(maindir, 'Demo_Code'))); %add all code to the path
%%
%load default parameters
[dirs, params] = getDefaultParameters(maindir);
load(fullfile(dirs.data2load, 'cell_metrics.mat'))
%create directory for saving figures
figdir = fullfile(maindir, 'Demo_Figures'); if ~isfolder(figdir); mkdir(figdir); end
%% use lap data, separate FAM and NOV
avg_trial = 1;
vel_fam = []; vel_nov = [];
lickrate_fam = []; lickrate_nov = [];
rate_fam = []; rate_nov = [];
stats_all = []; % for histogram, columns in order: R^2, F-stat, p-val, error variance. Each row is a neuron.
filelist = dir(fullfile(dirs.data2load, 'singlesess_ratemap_laps/*mat'));
for iFile = 1:length(filelist)
    animal = strsplit(filelist(iFile).name, '_');
    animal = animal{1}(2:end);
    if ~ismember(str2num(animal), params.WTmice)
        continue
    end
    load(fullfile(dirs.data2load, 'singlesess_ratemap_laps', filelist(iFile).name))
    sessstr = strsplit(filelist(iFile).name, '_'); sessstr = strcat(sessstr(1), {'_'}, sessstr(2)); sessstr = [sessstr{:} '_CA3'];
    celltype = cell_metrics.putativeCellType(find(ismember(cell_metrics.sessionName, sessstr)));
    ns_idx = find(ismember(celltype, 'Narrow Interneuron') & ~ismember(find(ismember(cell_metrics.sessionName, sessstr)), cell_metrics.tags.Bad));
    if isempty(ns_idx)
        continue
    end
    % get raw speed and firing rate. calculate stats for each neuron. 
    vel_fam_currday = outmap.bin5.fam.smoothSpeed; vel_nov_currday = outmap.bin5.nov.smoothSpeed;
    lickrate_fam_currday = outmap.bin5.fam.smoothLickrate; lickrate_nov_currday = outmap.bin5.nov.smoothLickrate;
    rate_fam_currday = outmap.bin5.fam.ratemap(:, :, ns_idx); rate_nov_currday = outmap.bin5.nov.ratemap(:, :, ns_idx);
    if avg_trial == 1
        % then avg across trials and neurons, concat, do normalization
        % and scaling for each day before plotting the entire data
        vel_fam_currday_avgtrial = nanmean(vel_fam_currday, 1); vel_nov_currday_avgtrial = nanmean(vel_nov_currday, 1);
        lickrate_fam_currday_avg_trial = nanmean(lickrate_fam_currday, 1); lickrate_nov_currday_avg_trial = nanmean(lickrate_nov_currday, 1);
        rate_fam_currday_avg_trialneuron = nanmean(rate_fam_currday, [1, 3]); rate_nov_currday_avg_trialneuron = nanmean(rate_nov_currday, [1, 3]);
        vel_fam = [vel_fam; vel_fam_currday_avgtrial];
        vel_nov = [vel_nov; vel_nov_currday_avgtrial];
        lickrate_fam = [lickrate_fam; lickrate_fam_currday_avg_trial];
        lickrate_nov = [lickrate_nov; lickrate_nov_currday_avg_trial];
        rate_fam = [rate_fam; rate_fam_currday_avg_trialneuron];
        rate_nov = [rate_nov; rate_nov_currday_avg_trialneuron];
    else % averaged spatial bins
        % then avg across spatial bins and neurons, concat, do normalization
        % and scaling for each day before plotting the entire data
        vel_fam_currday_avgbin = nanmean(vel_fam_currday, 2); vel_nov_currday_avgbin = nanmean(vel_nov_currday, 2);
        lickrate_fam_currday_avg_bin = nanmean(lickrate_fam_currday, 2); lickrate_nov_currday_avg_bin = nanmean(lickrate_nov_currday, 2);
        rate_fam_currday_avg_binneuron = nanmean(rate_fam_currday, [2, 3]); rate_nov_currday_avg_binneuron = nanmean(rate_nov_currday, [2, 3]);

        vel_fam = [vel_fam; vel_fam_currday_avgbin];
        vel_nov = [vel_nov; vel_nov_currday_avgbin];
        lickrate_fam = [lickrate_fam; lickrate_fam_currday_avg_bin];
        lickrate_nov = [lickrate_nov; lickrate_nov_currday_avg_bin];
        rate_fam = [rate_fam; rate_fam_currday_avg_binneuron];
        rate_nov = [rate_nov; rate_nov_currday_avg_binneuron];
    end
end
% plot
figure;
normalize = 1; scale = 1;
if normalize == 1
    rate_fam = (rate_fam - min(rate_fam, [], 2)) ./ (max(rate_fam, [], 2) - min(rate_fam, [], 2));
    rate_nov = (rate_nov - min(rate_nov, [], 2)) ./ (max(rate_nov, [], 2) - min(rate_nov, [], 2));
    vel_fam = (vel_fam - min(vel_fam, [], 2)) ./ (max(vel_fam, [], 2) - min(vel_fam, [], 2));
    vel_nov = (vel_nov - min(vel_nov, [], 2)) ./ (max(vel_nov, [], 2) - min(vel_nov, [], 2));
    lickrate_fam = (lickrate_fam - min(lickrate_fam, [], 2)) ./ (max(lickrate_fam, [], 2) - min(lickrate_fam, [], 2));
    lickrate_nov = (lickrate_nov - min(lickrate_nov, [], 2)) ./ (max(lickrate_nov, [], 2) - min(lickrate_nov, [], 2));
end
if scale == 1
    rate_fam = rate_fam - mean(rate_fam(:, 1:2), 2);
    rate_nov = rate_nov - mean(rate_nov(:, 1:2), 2);
    vel_fam = vel_fam - mean(vel_fam(:, 1:2), 2);
    vel_nov = vel_nov - mean(vel_nov(:, 1:2), 2);
    lickrate_fam = lickrate_fam - mean(lickrate_fam(:, 1:2), 2);
    lickrate_nov = lickrate_nov - mean(lickrate_nov(:, 1:2), 2);
end
% plot
t = tiledlayout(2, 2);
nexttile;
mdl = fitlm(reshape(vel_fam, [], 1), reshape(rate_fam, [], 1)); 
ph = plot(mdl);
set(ph(1), 'Color',params.colors_fam(end, :)); set(ph(2),'Color',params.colors_fam(end, :),'LineWidth',2); legend('off')
xlabel('norm velocity'); ylabel('norm NS FR');title(['r^2=' num2str(mdl.Rsquared.Adjusted) ',slope=' num2str(mdl.Coefficients.Estimate(2)), ',p=', num2str(mdl.ModelFitVsNullModel.Pvalue)]); ylim([-1 1]); xlim([-1 1])

nexttile;
mdl = fitlm(reshape(vel_nov, [], 1), reshape(rate_nov, [], 1)); 
ph = plot(mdl);
set(ph(1), 'Color',params.colors_nov(end-1, :)); set(ph(2), 'Color',params.colors_nov(end, :), 'LineWidth',2); legend('off')
xlabel('norm velocity'); ylabel('norm NS FR');title(['r^2=' num2str(mdl.Rsquared.Adjusted) ',slope=' num2str(mdl.Coefficients.Estimate(2)), ',p=', num2str(mdl.ModelFitVsNullModel.Pvalue)]); ylim([-1 1]); xlim([-1 1])

nexttile;
mdl = fitlm(reshape(lickrate_fam, [], 1), reshape(rate_fam, [], 1));
ph = plot(mdl);
set(ph(1), 'Color',params.colors_fam(end, :)); set(ph(2), 'Color',params.colors_fam(end, :), 'LineWidth',2); legend('off')
xlabel('norm lickrate'); ylabel('norm NS FR');title(['r^2=' num2str(mdl.Rsquared.Adjusted) ',slope=' num2str(mdl.Coefficients.Estimate(2)), ',p=', num2str(mdl.ModelFitVsNullModel.Pvalue)]); ylim([-1 1]); xlim([-1 1])

nexttile;
mdl = fitlm(reshape(lickrate_nov, [], 1), reshape(rate_nov, [], 1));
ph = plot(mdl);
set(ph(1), 'Color',params.colors_nov(end-1, :)); set(ph(2), 'Color',params.colors_nov(end, :), 'LineWidth',2); legend('off')
xlabel('norm lickrate'); ylabel('norm NS FR');title(['r^2=' num2str(mdl.Rsquared.Adjusted) ',slope=' num2str(mdl.Coefficients.Estimate(2)), ',p=', num2str(mdl.ModelFitVsNullModel.Pvalue)]); ylim([-1 1]); xlim([-1 1])

makefigurepretty(gcf,1)
figname = 'ExtendedDataFigure02_FG';
savefigALP([figdir '/'], figname, 'filetype', 'pdf')

%% compare linear with nonlinear methods and print out result
% Number of folds for cross-validation
behavdata = {vel_fam, vel_nov, lickrate_fam, lickrate_nov};
neuraldata = {rate_fam, rate_nov, rate_fam, rate_nov};
% Preallocate arrays for R-squares
k = 5; 
r2_linear = zeros(k, 4);
r2_rf = zeros(k, 4);
r2_svm = zeros(k, 4);
r2_knn = zeros(k, 4);
r2_gpr = zeros(k, 4);
for iBehavType = 1:4
    X = behavdata{iBehavType}; X = reshape(X, [], 1);
    Y = neuraldata{iBehavType}; Y = reshape(Y, [], 1);
    
    cv = cvpartition(length(X), 'KFold', k);
    
    for i = 1:k
        % Get training and test indices
        train_idx = training(cv, i);
        test_idx = test(cv, i);
    
        % Split data into training and test sets
        X_train = X(train_idx);
        Y_train = Y(train_idx);
        X_test = X(test_idx);
        Y_test = Y(test_idx);
        
        % Linear model
        lm = fitlm(X_train, Y_train);
        Y_pred_lm = predict(lm, X_test);
        r2_linear(i, iBehavType) = 1 - sum((Y_test - Y_pred_lm).^2) / sum((Y_test - mean(Y_test)).^2);
        
        % Random forest regression
        rf_model = TreeBagger(100, X_train, Y_train, 'Method', 'regression');
        Y_pred_rf = predict(rf_model, X_test);
        % Y_pred_rf = cell2mat(Y_pred_rf); % Convert to numeric array
        r2_rf(i, iBehavType) = 1 - sum((Y_test - Y_pred_rf).^2) / sum((Y_test - mean(Y_test)).^2);
    
        % Support Vector Machine regression
        svr_model = fitrsvm(X_train, Y_train, 'KernelFunction', 'gaussian');
        Y_pred_svm = predict(svr_model, X_test);
        r2_svm(i, iBehavType) = 1 - sum((Y_test - Y_pred_svm).^2) / sum((Y_test - mean(Y_test)).^2);

        % K-Nearest Neighbors regression
        knn_model = fitcknn(X_train, Y_train, 'NumNeighbors', 5, 'Distance', 'euclidean');
        Y_pred_knn = predict(knn_model, X_test);
        r2_knn(i, iBehavType) = 1 - sum((Y_test - Y_pred_knn).^2) / sum((Y_test - mean(Y_test)).^2);

        % Gaussion Process regression
        gpr_model = fitrgp(X_train, Y_train);
        Y_pred_gpr = predict(gpr_model, X_test);
        r2_gpr(i, iBehavType) = 1 - sum((Y_test - Y_pred_gpr).^2) / sum((Y_test - mean(Y_test)).^2);

    end
end

% Average R-squares
r2_linear_mean = mean(r2_linear);
r2_rf_mean = mean(r2_rf);
r2_svm_mean = mean(r2_svm);
r2_knn_mean = mean(r2_knn);
r2_gpr_mean = mean(r2_gpr);

fprintf('Cross-validated R-square for linear model: %.4f\n', r2_linear_mean);
fprintf('Cross-validated R-square for random forest: %.4f\n', r2_rf_mean);
fprintf('Cross-validated R-square for SVM regression: %.4f\n', r2_svm_mean);
fprintf('Cross-validated R-square for KNN regression: %.4f\n', r2_knn_mean);
fprintf('Cross-validated R-square for Gaussion Process regression: %.4f\n', r2_gpr_mean);
disp('comparisons: velocity_fam, velocity_nov, lickrate_fam, lickrate_nov')


