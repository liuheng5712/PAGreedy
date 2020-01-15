function [selected_UNI, time_UNI, record_indices] = UNI(index, k, b, beta, gamma)

%input: index:the dataset, k is amount of features need to be selected, b is the bootstrapping indices
%beta and gamma should be between zero and one
%output: seleced features and time matrices
% index = 19;
% k = 100;
% beta = 0.5;
% gamma = 0.5;

fid = fopen('/home/hengl/matlab/bin/scripts/DAJMI/names2');
file_name = {};
for i=1:30
    line = fgets(fid);
    if ischar(line)
        file_name{i} = line;
    end
end
expression = '[A-Za-z0-9-_]*';
for i=1:30
    str = regexp(file_name{i}, expression);
    file_name{i} = file_name{i}(1: str(1,2)-2);
end
path = strcat('/home/hengl/matlab/bin/scripts/DAJMI/data/', file_name{index}, '.mat');
file = load(path);
data = struct2cell(file);
data = data{1};
f = size(data,2)-1;
data = data(b,:);
y = data(:, size(data,2));

reles = [];   %keeep the mutual information terms w.r.t. label
selected_UNI = []; %keep the selected feature subset
current = 0; %keep track of the current selected feature 
al = 30;
record_indices = [];
record_score = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ma = 0;
for i=1:f
    re = mi(data(:,i),y);
    reles = [reles re];
    if ma < re
        ma = re;
        current = i;
    end
end
feature_set = 1:f; %To get rid of the all-zero columns

selected_UNI = [selected_UNI current];
feature_set = feature_set(feature_set ~= current);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_UNI = [0];

start = tic;
for q=1:k-1  %select k features
    score = [];
    
    ma = 0;
    pin = 0;     %to keep track of the feature that is about to select
    for i=feature_set
        temp = 0;    %to keep track of the current sum of culmuative joint mutual information
        for j=selected_UNI
            temp = temp + beta*mi(data(:,i), data(:,j)) - gamma*cmi(data(:,i), data(:,j), y);
        end
        temp = reles(i) - temp;
        score = [score temp];
        
        if ma < temp
            ma = temp;
            pin = i;
        end
    end
    
    
    %find the the diagonal pattern
    [B I] = sort(score, 'descend');
    score = B;
    order = feature_set(I);
    record_indices(q,:) = order(1:al);
    record_score(q,:) = score(1:al);
    
    %selecting
    selected_UNI = [selected_UNI pin];
    feature_set = feature_set(feature_set ~= pin);
    time_UNI = [time_UNI toc(start)];
    disp('selecting...')
    disp(size(selected_UNI,2))
end     

