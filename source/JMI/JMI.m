function [selected_JMI, time_JMI, record_indices] = JMI(index, k, b)

%input: index:the dataset, k is amount of features need to be selected, b is the bootstrapping indices
%output: seleced features and time matrices
% index = 21;
% k = 100;

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


selected_JMI = []; %keep the selected feature subset
current = 0; %keep track of the current selected feature 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ma = 0;
for i=1:f
    re = mi(data(:,i),y);
    if ma < re
        ma = re;
        current = i;
    end
end
feature_set = [1:f]; %To get rid of the all-zero columns

previous_selected = current; %use this to keep track of previous selected feature
selected_JMI = [selected_JMI current];
feature_set = feature_set(feature_set ~= current);
score = zeros(1,f-1); %store the score of a candidate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ma = 0; %keep track of the max of the culmuative joint mutual information
time_JMI = [0];

al = 30;
record_indices = [];
record_score = [];

start = tic;
for q=1:k-1  %select k features

    score = [];
    ma = 0;
    pin = 0;     %to keep track of the feature that is about to select
    for i=feature_set
        temp = 0;

        for j = selected_JMI
            jo = joint([data(:,i),data(:,j)],[]);
            temp = temp + mi(jo,y);
        end
        score = [score temp];
        if ma < temp
            ma = temp;
            pin = i;
        end
    end
    
    %find the the diagonal pattern
    [B,I] = sort(score, 'descend');
    credits = B;
    order = feature_set(I);
    record_indices(q,:) = order(1:al);
    record_score(q,:) = credits(1:al);
    
    
    selected_JMI = [selected_JMI pin];
    score = score(feature_set ~= pin);
    feature_set = feature_set(feature_set ~= pin);
    previous_selected = pin;
    time_JMI = [time_JMI toc(start)];
    
    %selecting
    disp('selecting...')
    disp(size(selected_JMI,2))
end     


