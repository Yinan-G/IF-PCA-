%%%% Load  Single-Cell data  
%Name set: (load p-by-n data matrix)
%camp1, 13111-by-777
%camp2, 11233-by-734
%darmanis, 13400-by-466
%deng, 16347-by-268
%goolam, 21199-by-124
%grun, 5547-by-1502
%li, 25369-by-561
%patel, 5948-by-430
name = 'grun';
path = "Dataset/sc_data"; 


filePath = fullfile(path, [name, '-x-filter.txt']);
SC_raw_data = dlmread(filePath, ',', 0, 0);
SC_raw_data = SC_raw_data';
[n, p] = size(SC_raw_data);

%SC_raw_data = bsxfun(@minus, SC_raw_data, min(SC_raw_data, [], 1)); % This step is for Patel dataset

file_y_Path = fullfile(path, [name, '-y.txt'])
if ~isreal(SC_raw_data)
    error('Input data must be a real array for built-in distances.');
end

% Load true class labels
labels_true = dlmread(file_y_Path, ' ', 0, 0);

% Calculate the number of classes
numClasses = max(labels_true) - min(labels_true) + 1;


knns=15;
knnf=50;
if p>10000
    knnf=100;
end

% %%%%%%%%%%% Load  Microarray data  
% % Name set: (load p-by-n data matrix)
% % brain, 5597-by-42
% % breast, 22215-by-276
% % colon, 2000-by-62
% % leukemia, 3571-by-72
% % lung1, 12533-by-181
% % lung2, 12600-by-203
% % lymphoma, 4026-by-62
% % prostate, 6033-by-102
% % srbct, 2308-by-63
% % su, 7909-by-174
% 
% name = 'colon';
% path = "Dataset/microarray_data"; 
% 
% filePath = fullfile(path, [name, '.x.txt']);
% 
% 
% SC_raw_data = dlmread(filePath, ',', 0, 0);
% SC_raw_data = SC_raw_data';
% file_y_Path = fullfile(path, [name, '.y.txt'])
% if ~isreal(SC_raw_data)
%     error('Input data must be a real array for built-in distances.');
% end
% 
% % Load true class labels
% labels_true = dlmread(file_y_Path, ' ', 0, 0);
% 
% % Calculate the number of classes
% numClasses = max(labels_true) - min(labels_true) + 1;
% 
% knns=10;
% knnf=10;



%%%%%%%%%%%

tic
% Transform data then Manifold denoise with denoising step in scAMF

fea_raw = full(SC_raw_data);
fea_raw_fit = manfit(fea_raw, knns);
fea_raw_fit = double(fea_raw_fit);
fea_raw_duo_fit = manfit(fea_raw_fit', knnf);
fea_raw_duo_fit = double(fea_raw_duo_fit');



fea_log =  transform(SC_raw_data,'log');
fea_log_fit = manfit(fea_log, knns);
fea_log_fit = double(fea_log_fit);
fea_log_duo_fit = manfit(fea_log_fit', knnf);
fea_log_duo_fit = double(fea_log_duo_fit');


%Save variables to MAT-file
f_name2 = '_Y2_denoised_data.mat';
f_name = [name, f_name2];

save(f_name, 'fea_raw_fit', 'fea_log_fit', 'fea_raw_duo_fit', 'fea_log_duo_fit');
%save(f_name, 'fea_raw_fit', 'fea_raw_duo_fit');

runtime = toc;

disp(['Runtime: ', num2str(runtime), ' seconds']);

