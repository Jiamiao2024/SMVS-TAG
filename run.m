clear;
clc;
warning off;
addpath(genpath('./'));

%% Dataset settings
ds = {
    'Dermatology', ...
    'Scene15', ...
    };

dsPath = fullfile(pwd, 'datasets');
resPath = fullfile(pwd, 'res');

if ~exist(resPath, 'dir')
    mkdir(resPath);
end

%% Run experiments
for dsi = 1:length(ds)

    dataName = ds{dsi};
    fprintf('\nDataset: %s\n', dataName);

    dataFile = fullfile(dsPath, [dataName, '.mat']);

    if ~exist(dataFile, 'file')
        error('Dataset file not found: %s', dataFile);
    end

    load(dataFile);

    k = length(unique(Y));

    %% Hyperparameter settings
    lambda_list = [1,10,100,1000];
    anchor_list = [1, 3, 5] * k;
    p_list = [0.1, 0.3, 0.5];
    beta_list = [1e-4, 1e-2, 1e-1, 1];

    %% Output file
    txtFile = fullfile(resPath, [dataName, '_result.txt']);
    fid = fopen(txtFile, 'w');
    if fid == -1
        error('Cannot create result file: %s', txtFile);
    end

    %% Write header
    fprintf(fid,'lambda\tanchor\tp\tbeta\tACC\tNMI\tPurity\tFscore\tPrecision\tRecall\tAR\tEntropy\tTime(s)\n');

    %% Parameter search
    for ll = 1:length(lambda_list)
        lambda = lambda_list(ll);
        for aa = 1:length(anchor_list)
            anchor = anchor_list(aa);
            for pp = 1:length(p_list)
                p = p_list(pp);
                for bb = 1:length(beta_list)
                    beta = beta_list(bb);
                    fprintf('Running: lambda=%.15g, anchor=%d, p=%.15g, beta=%.15g\n', lambda, anchor, p, beta);
                    tic;
                    [U, V, VW, W, Z, alpha, iter, obj, Result] = algo_chd(X, Y, lambda, anchor, k, p, beta);
                    res = myNMIACCwithmean(U, Y, k);
                    timer = toc;

                    %% Convert result to row vector
                    res = res(:)';
                    
                    %% Check metric length
                    if length(res) ~= 8
                        error('The evaluation result must contain 8 metrics.');
                    end

                    %% Write results
                    fprintf(fid, '%.15g\t%d\t%.15g\t%.15g', lambda, anchor, p, beta);
                    fprintf(fid, ...
                        '\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f', ...
                        res(1), res(2), res(3), res(4), res(5), res(6), res(7), res(8));
                    fprintf(fid, '\t%.4f\n', timer);
                    fprintf('ACC=%.4f, NMI=%.4f, Purity=%.4f, Fscore=%.4f, Time=%.4f s\n', ...
                        res(1), res(2), res(3), res(4), timer);

                end
            end
        end
    end

    fclose(fid);
    fprintf('Results saved to: %s\n', txtFile);
    clear X Y;
end