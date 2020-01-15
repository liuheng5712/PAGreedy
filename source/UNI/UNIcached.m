%function [selected_UNIcached, time_UNIcached, record_indices] = UNIcached(index, k, b)

%input: index:the dataset, k is amount of features need to be selected, b is the bootstrapping indices
%output: seleced features and time matrices
index = 9;
k = 101;
beta = 0.5;
gamma = 0.5;
b = 1:452;

fid = fopen('/home/hengl/matlab/bin/scripts/PAGreedy2.0/names2');
file_name = {};
for i=1:10
    line = fgets(fid);
    if ischar(line)
        file_name{i} = line;
    end
end
expression = '[A-Za-z0-9-_]*';
for i=1:10
    str = regexp(file_name{i}, expression);
    file_name{i} = file_name{i}(1: str(1,2)-2);
end
path = strcat('/home/hengl/matlab/bin/scripts/PAGreedy2.0/data/', file_name{index}, '.mat');
file = load(path);
data = struct2cell(file);
data = data{1};
f = size(data,2)-1;
data = data(b,:);
y = data(:, size(data,2));

reles = [];   %keeep the mutual information terms w.r.t. label
selected_UNIcached = []; %keep the selected feature subset
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
feature_set = 1:f;

selected_UNIcached = [selected_UNIcached current];
score_vector = reles(feature_set ~= current);
feature_set = feature_set(feature_set ~= current);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_UNIcached = [0];
start = tic;
previous = current; %keep track of previously selected feature
for q=1:k-1  %select k features
    
    ma = -10000;
    pin = 0;     %to keep track of the feature that is about to select
    for i = feature_set

        score_vector(feature_set == i) =  score_vector(feature_set == i) - beta*mi(data(:,i),data(:,previous)) + gamma*cmi(data(:,i), data(:,previous), y);
        if ma < score_vector(feature_set == i)
            ma = score_vector(feature_set == i);
            pin = i;
        end
    end
    
    
    %find the the diagonal pattern
    [B I] = sort(score_vector, 'descend');
    score = B;
    order = feature_set(I);
    record_indices(q,:) = order(1:al);
    record_score(q,:) = score(1:al);
    
    %selecting
    selected_UNIcached = [selected_UNIcached pin];
    previous = pin;
    score_vector = score_vector(feature_set ~= pin);
    feature_set = feature_set(feature_set ~= pin);
    time_UNIcached = [time_UNIcached toc(start)];
    disp('selecting...')
    disp(size(selected_UNIcached,2))
end     

