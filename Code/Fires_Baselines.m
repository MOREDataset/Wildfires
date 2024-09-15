%% TASK1 BASELINES: Season-based naive model
clear
clc
fprintf('Task1: Season-based naive baseline\nHrz;Accuracy;FAR\n');
for i = [0 7 30 90]
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireDetection_trainsetV4_' num2str(i) 'Days.csv']);
    opts = detectImportOptions(path1);
    trainset = readtable(path1,opts);

    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireDetection_testsetV4_' num2str(i) 'Days.csv']);
    opts = detectImportOptions(path1);
    testset = readtable(path1,opts);

    % Find most frequent season only for fire events
    trainset_posClass = trainset(trainset.isFire==1,:);
    Task1_NaiveSeason = mode((trainset_posClass.S));

    % Define targets and predictions
    y_test = testset.isFire;

    y_preds = zeros(size(y_test));
    [rows,~] = find(testset.S==Task1_NaiveSeason);
    y_preds(rows) = ones(size(rows));

    % Compute performance
    [C,~] = confusionmat(y_test, y_preds);
    TP = C(1,1);
    TN = C(2,2);
    FP = C(2,1);
    FN = C(1,2);

    accuracy = (TP + TN) / (TP + TN + FP + FN);
    FAR = FP / (FP + TN);

    % Display results
    fprintf('-;-;-;-;-;%.3f%%; %.3f%%;', accuracy * 100, FAR* 100);
    % fprintf('%d; %.3f%%; %.3f%%\n', i, accuracy * 100, FAR* 100);
end

%% TASK1 BASELINES: Location-based naive model
clc
fprintf('Task1: Top 10%% Frequent Locations-based naive baseline\nHrz;Accuracy;FAR\n');
for i = [0 7 30 90]
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireDetection_trainsetV4_' num2str(i) 'Days.csv']);
    opts = detectImportOptions(path1);
    trainset = readtable(path1,opts);

    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireDetection_testsetV4_' num2str(i) 'Days.csv']);
    opts = detectImportOptions(path1);
    testset = readtable(path1,opts);

    % Find top 10% frequent location (avgLatitude, avgLongitude pairs) only for fire events
    trainset_posClass = trainset(trainset.isFire == 1, :);
    locPairs = strcat(string(trainset_posClass.avgLatitude), '_', string(trainset_posClass.avgLongitude));
    [G,locPairsGroup] = findgroups(locPairs);
    % locCounts = splitapply(@numel, locPairs, G);
    % [sortedLocCounts, sortIdx] = sort(locCounts, 'descend');
    % topLocPairs = locPairsGroup(sortIdx);
    % 
    % % Calculate the top 10% index
    % top10PercentIdx = ceil(numel(topLocPairs) * 0.2);
    % topLocPairs = topLocPairs(1:top10PercentIdx);

    % Define targets and predictions
    y_test = testset.isFire;
    y_preds = zeros(size(y_test));

    for j = 1:height(testset)
        % Check if the location is in the top 10% frequent locations
        locPair = strcat(string(testset.avgLatitude(j)), '_', string(testset.avgLongitude(j)));
        isLocation = any(~ismember(locPairsGroup, locPair));%topLocPairs
        if isLocation == 1
            y_preds(j) = 1;
        end
    end

    % Compute performance
    [C, ~] = confusionmat(y_test, y_preds);
    TP = C(1,1);
    TN = C(2,2);
    FP = C(2,1);
    FN = C(1,2);

    accuracy = (TP + TN) / (TP + TN + FP + FN);
    FAR = FP / (FP + TN);

    % Display results
    fprintf('-;-;-;-;-;%.3f%%; %.3f%%;', accuracy * 100, FAR* 100);
    % fprintf('%d; %.3f%%; %.3f%%\n', i, accuracy * 100, FAR * 100);
end

%% TASK1 BASELINES: Blended Season and Location-based naive model
clc
fprintf('Task1: Blended Season and Location-based naive baseline\nHrz;Accuracy;FAR\n');
for i = [0 7 30 90]
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireDetection_trainsetV4_' num2str(i) 'Days.csv']);
    opts = detectImportOptions(path1);
    trainset = readtable(path1,opts);

    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireDetection_testsetV4_' num2str(i) 'Days.csv']);
    opts = detectImportOptions(path1);
    testset = readtable(path1,opts);

    % Season-based Approach
    trainset_posClass_season = trainset(trainset.isFire == 1, :);
    mostFrequentSeason = mode(trainset_posClass_season.S);

    % Location-based Approach
    trainset_posClass_location = trainset(trainset.isFire == 1, :);
    locPairs = strcat(string(trainset_posClass_location.avgLatitude), '_', string(trainset_posClass_location.avgLongitude));
    [G,locPairsGroup] = findgroups(locPairs);
    % locCounts = splitapply(@numel, locPairs, G);
    % [sortedLocCounts, sortIdx] = sort(locCounts, 'descend');
    % topLocPairs = locPairsGroup(sortIdx);
    % top10PercentIdx = ceil(numel(topLocPairs) * 0.1);
    % topLocPairs = topLocPairs(1:top10PercentIdx);

    % Define targets and predictions
    y_test = testset.isFire;
    y_preds = zeros(size(y_test));

    for j = 1:height(testset)
        % Check if the location is in the top 10% frequent locations
        locPair = strcat(string(testset.avgLatitude(j)), '_', string(testset.avgLongitude(j)));
        isLocation = any(~ismember(locPairsGroup, locPair));%topLocPairs

        % Check if the season matches the most frequent fire season
        isTopSeason = testset.S(j) == mostFrequentSeason;

        % Predict fire if either condition is true
        if isLocation || isTopSeason
            y_preds(j) = 1;
        end
    end

    % Compute performance
    [C, ~] = confusionmat(y_test, y_preds);
    TP = C(1,1);
    TN = C(2,2);
    FP = C(2,1);
    FN = C(1,2);

    accuracy = (TP + TN) / (TP + TN + FP + FN);
    FAR = FP / (FP + TN);

    % Display results
    fprintf('%d; %.3f%%; %.3f%%\n', i, accuracy * 100, FAR * 100);
end
%% TASK2 BASELINES: Season-based naive model
clear 
clc
fprintf('Task2: Season-based naive baseline\nHrz;Accuracy;FAR\n');
for i = [0 7 30 90]
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireBurntArea_trainset_' num2str(i) 'Days.csv']);
    opts = detectImportOptions(path1);
    trainset = readtable(path1,opts);

    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireBurntArea_testset_' num2str(i) 'Days.csv']);
    opts = detectImportOptions(path1);
    testset = readtable(path1,opts);

    % Find most frequent season
    [freq_S0,~] = find((trainset.S==0));
    [freq_S1,~] = find((trainset.S==1));
    [freq_S2,~] = find((trainset.S==2));
    [freq_S3,~] = find((trainset.S==3));
    Task3_NaiveSeason0 = mean(trainset.totalArea_km2(freq_S0));
    Task3_NaiveSeason1 = mean(trainset.totalArea_km2(freq_S1));
    Task3_NaiveSeason2 = mean(trainset.totalArea_km2(freq_S2));
    Task3_NaiveSeason3 = mean(trainset.totalArea_km2(freq_S3));
    
    % Define targets and predictions
    y_train = trainset.totalArea_km2;
    y_test = testset.totalArea_km2;

    y_preds = zeros(size(y_test));
    [rows_S0,~] = find(testset.S==0);
    [rows_S1,~] = find(testset.S==1);
    [rows_S2,~] = find(testset.S==2);
    [rows_S3,~] = find(testset.S==3);
    y_preds(rows_S0) = Task3_NaiveSeason0*ones(size(rows_S0));
    y_preds(rows_S1) = Task3_NaiveSeason1*ones(size(rows_S1));
    y_preds(rows_S2) = Task3_NaiveSeason2*ones(size(rows_S2));
    y_preds(rows_S3) = Task3_NaiveSeason3*ones(size(rows_S3));

    % Compute performance
    rmse = sqrt(mean((y_test-y_preds).^2));
    mase = mean_absolute_scaled_error(y_test, y_preds, y_train);

    % Display results
    fprintf('-;-;-;-;-;%.3f; %.3f;', rmse, mase);
    % fprintf('%d; %.3f; %.3f\n', i, rmse, mase);
end

%% TASK3 BASELINES: Season-based naive model
clear 
clc
fprintf('Task3: Season-based naive baseline\nHrz;Accuracy;FAR\n');
for i = [0 7 30 90]
    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireDuration_trainset_' num2str(i) 'Days.csv']);
    opts = detectImportOptions(path1);
    trainset = readtable(path1,opts);

    path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\Fire\NewWork\V2',['FireDuration_testset_' num2str(i) 'Days.csv']);
    opts = detectImportOptions(path1);
    testset = readtable(path1,opts);

    % Find most frequent season
    [freq_S0,~] = find((trainset.S==0));
    [freq_S1,~] = find((trainset.S==1));
    [freq_S2,~] = find((trainset.S==2));
    [freq_S3,~] = find((trainset.S==3));
    Task3_NaiveSeason0 = mean(trainset.Duration(freq_S0));
    Task3_NaiveSeason1 = mean(trainset.Duration(freq_S1));
    Task3_NaiveSeason2 = mean(trainset.Duration(freq_S2));
    Task3_NaiveSeason3 = mean(trainset.Duration(freq_S3));
    
    % Define targets and predictions
    y_train = trainset.Duration;
    y_test = testset.Duration;

    y_preds = zeros(size(y_test));
    [rows_S0,~] = find(testset.S==0);
    [rows_S1,~] = find(testset.S==1);
    [rows_S2,~] = find(testset.S==2);
    [rows_S3,~] = find(testset.S==3);
    y_preds(rows_S0) = Task3_NaiveSeason0*ones(size(rows_S0));
    y_preds(rows_S1) = Task3_NaiveSeason1*ones(size(rows_S1));
    y_preds(rows_S2) = Task3_NaiveSeason2*ones(size(rows_S2));
    y_preds(rows_S3) = Task3_NaiveSeason3*ones(size(rows_S3));

    % Compute performance
    rmse = sqrt(mean((y_test-y_preds).^2));
    mase = mean_absolute_scaled_error(y_test, y_preds, y_train);

    % Display results
    fprintf('-;-;-;-;-;%.3f; %.3f;', rmse, mase);
    % fprintf('%d; %.3f; %.3f\n', i, rmse, mase);
end




%% functions

function mae = mean_absolute_error(y_true, y_pred)
    % Calculate the mean absolute error
    mae = mean(abs(y_true - y_pred));
end


function mase = mean_absolute_scaled_error(y_true, y_pred, y_train)
    % Compute the Mean Absolute Scaled Error (MASE)

    % Check for the length of y_train
    if length(y_train) < 2
        error('Length of y_train should be at least 2.');
    end

    % Calculate MAE for the predictions
    mae = mean_absolute_error(y_true, y_pred);

    % Calculate the MAE of the one-step naive forecast method
    naive_forecast = y_train(1:end-1);
    naive_true = y_train(2:end);
    mae_naive = mean_absolute_error(naive_true, naive_forecast);

    % Handle the case when naive MAE is zero to avoid division by zero
    if mae_naive == 0
        mase = inf;
        if mae == 0
            mase = 0;
        end
        return;
    end

    % Calculate MASE
    mase = mae / mae_naive;
end
