%function [selected_PAJMImerged, time_PAJMImerged, window_trend] = PAJMIcachedmerged(dataindex, k, b, M)

%input: dataindex is the data, k is the amount of features selected, b is the bootstrapping indices, M is amount of worker available
%output: selected features and time vector
dataindex = 10;
k = 200;
M = 10;
b = 1:6238;

%Data preparation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
path = strcat('/home/hengl/matlab/bin/scripts/PAGreedy2.0/data/', file_name{dataindex}, '.mat');
file = load(path);
data = struct2cell(file);
data = data{1};
f = size(data,2)-1;
data = data(b,:);
y = data(:, size(data,2));


%Initialization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_PAJMImerged = [0];
selected_PAJMImerged = []; %keep the selected feature subset
current = 0; %keep track of the current selected feature 
ma = 0;

for i=1:f
    re = mi(data(:,i),y);
    if ma < re
        ma = re;
        current = i;
    end
end
feature_set = [1:f]; %To get rid of the all-zero columns
selected_PAJMImerged = [selected_PAJMImerged current];
feature_set = feature_set(feature_set ~= current);

parpool(4); %start the parallel pool and distributing the feature space
poolobj = gcp;
addAttachedFiles(poolobj,{'mi.m', 'joint.m'});


%Prepare for parallelzation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dist_feature_set = [];   %the rows is the indices for worker
sizeBlock = ceil(size(feature_set,2)/M);
parfor i = 1:M-1
    index = ((i-1)*sizeBlock+1):i*sizeBlock;
    dist_feature_set(i,:) = feature_set(index);
end
index = ((M-1)*sizeBlock+1):size(feature_set,2);
dist_feature_set(M,:) = [feature_set(index) zeros(1,(sizeBlock - size(index,2)))]; %denote the empty spots as 0
score_dist_features = zeros(M, sizeBlock); %this is used to store the scores for each candidate
previous_selected = selected_PAJMImerged; %this keeps track of selected subset previous step

start = tic; %Start the clock

window_trend = [];
selected_trend = [];
%selecting k features%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while size(selected_PAJMImerged,2) <= k %select k or more features  
    %calculating the score of each candidate feature
    parfor i = 1:M  
        features = dist_feature_set(i,:);
        scores = score_dist_features(i,:); %old scores
        temp_scores = [];
        for p = features
            temp = 0;
            if p ~= 0 %we denote the selected/unavailable features as 0
                for q = previous_selected
                   temp = temp + mi(joint([data(:,p),data(:,q)],[]),y);
                end
                temp_scores = [temp_scores temp];
            else
                temp_scores = [temp_scores -1];
            end
        end
        score_dist_features(i,:) = scores + temp_scores;
    end
    
    %sort the candidate features and get the sequence
    scores = score_dist_features'; 
    scores = scores(:);
    [B, I] = sort(scores, 'descend');
    sorted_features = feature_set(I(B >= 0));
    B = B(B >= 0);

    %for each candidate calculate its boundary and verify the candidates in parallel 
    flag = ones(1, M+1); %1 means no intrusion
    intruders = zeros(1, M+1); %mark the intruders
    
    MIcorrespond = []; %prepare to store the temperory MI's 
    cell = {};                %the cached comparisons
    windows = [];
    summation = [];
    
    parfor i = 2:M+1   %find the interfering zone
        boundary = 0;
        pre = []; %joint informativeness with features located ahead
        
        for s = 1:i-1 %find the boundary of the interfering zone 
            temp = mi(joint([data(:,sorted_features(i)),data(:,sorted_features(s))],[]), y);
            pre = [pre temp];
            boundary = boundary + temp - h(y); 
        end
        
        summation(i) = sum(pre);
        windowsize = 1;
        for p = i+1:size(sorted_features,2)
            if B(i) - B(p) + boundary < 0
                windowsize = windowsize + 1;
            else
                break;
            end
        end   
        windowsize = windowsize + i;
        if windowsize  > size(sorted_features,2)
            windowsize = size(sorted_features,2)
        end
        windows = [windows windowsize];
        
        cell{i} = B(i) + sum(pre); 
    end
    window_trend = [window_trend max(windows)];
    windows = max(windows);
    
    
    %here we add this procedure to calculate all the information theoretic
    %terms to avoid the redundant calculation of previous implementation
    parfor i = 1:M
        for j = i+1:windows
            temp = mi(joint([data(:,sorted_features(i)),data(:,sorted_features(j))],[]), y);
            MIcorrespond = [MIcorrespond; [sorted_features(i) sorted_features(j) temp]];
        end
    end
    MImatrix = sparse([MIcorrespond(:,1);MIcorrespond(:,2)], [MIcorrespond(:,2);MIcorrespond(:,1)],[MIcorrespond(:,3);MIcorrespond(:,3)]); %store the useful MI's and make it symmetric
    
    
    %here we do the identification
    parfor i = 2:M+1   
        temperory = [];
        intrude_max_objective = 0; % to keep track of the score of the strongest intruder
        for p = i+1:windows        %iterate all possible intruders until the boundary reaches
            objective = B(i) - B(p) + summation(i);
            right = B(p);
            for q = 1:i-1  %each candidate is evaluated w.r.t. its predecessor 
                temp = MImatrix(sorted_features(p), sorted_features(q));
                right = right + temp;
                objective =  objective - temp;
            end
            temperory = [temperory; [cell{i} right]];
            if objective < 0 
                flag(i) = 0;
                if (B(i) + summation(i) - objective) > intrude_max_objective;
                    intruders(i) = sorted_features(p);
                    intrude_max_objective = B(i) + summation(i) - objective;
                end
            end 
        end 
        cell{i} = temperory;
    end
    disp(sorted_features(1:M+1))  %**testing**%
    disp(flag) %**testing**%
    disp(intruders)  %**testing**%
        
    
    
    %merging the vectors to select candidates
    previous_selected = [];
    for i = 1:size(flag,2) %merge and shuffle the flag vector
        if flag(i) == 1 %it's a Hit
            previous_selected = [previous_selected sorted_features(i)];
        else  %it's a Miss
            if size(find(intruders(i) == sorted_features(1:M)),2) > 0 %the intruder is within the window
                newplace = find(sorted_features == intruders(i));
                originalfront = sorted_features(i); originalrear = sorted_features(newplace);
                difference = size(cell{i},1) - size(cell{newplace},1); oldcell = cell{i}; cell{newplace} = oldcell(difference+1:size(cell{i},1),:);
                sorted_features(sorted_features == intruders(i)) = sorted_features(i); sorted_features(i) = intruders(i);%switching
                flag(i) = 1; flag(i+1:newplace) = 0; intruders(i+1:newplace) = 0;
                previous_selected = [previous_selected sorted_features(i)];
                disp(sorted_features(1:M+1))  %**testing**%
                for t = i+1:newplace  %updating comparisons
                    comparison = cell{t};
                    if t < newplace
                        comparison(:,1) = comparison(:,1) - MImatrix(originalfront, sorted_features(t)) + MImatrix(originalrear, sorted_features(t));
                        parfor z = 1:size(comparison,1)
                            if (z + t) == newplace
                                tempcell = cell{i};
                                comparison(z,2) = tempcell(1,1);
                                for w = i:t-1
                                    comparison(z,2) = comparison(z,2) + MImatrix(sorted_features(w), originalfront);
                                end
                            else
                                comparison(z,2) = comparison(z,2) - MImatrix(originalfront, sorted_features(t+z)) + MImatrix(originalrear, sorted_features(t+z));
                            end
                        end
                    else
                        for m = i:newplace-1
                            for z = 1:size(comparison,1)
                                comparison(z,1) = comparison(z,1) + MImatrix(originalfront, sorted_features(m));
                                comparison(z,2) = comparison(z,2) + MImatrix(sorted_features(newplace+z), sorted_features(m));
                            end
                        end
                    end
                    cell{t} = comparison;
                    temp = find(comparison(:,1) <= comparison(:,2));
                    if size(temp,1) == 0
                        flag(t) = 1;
                    else
                        flag(t) = 0;
                        tempintruder = sorted_features(t + find(comparison(:,2) == max(comparison(:,2))));
                        intruders(t) = tempintruder(1);
                    end
                end
                disp(flag)  %**testing**%
                disp(intruders)  %**testing**%
            else
                break;
            end
        end       
    end
    
    
    index = find(flag == 0);
    if size(index,2) > 0 %absorb the Maximal Intruder
        previous_selected = [previous_selected intruders(index(1))];
    end

    parfor i = 1:M
        features = dist_feature_set(i,:);
        scores = score_dist_features(i,:);
        index = 0;
        for j = previous_selected
            index = find(features == j);
            features(index) = 0;
            scores(index) = -1;
        end
        dist_feature_set(i,:) = features;
        score_dist_features(i,:) = scores;
    end
    selected_PAJMImerged = [selected_PAJMImerged previous_selected];
    selected_trend = [selected_trend size(selected_PAJMImerged,2)];
    
    disp('selecting...')
    disp(size(selected_PAJMImerged,2))

    
    time_PAJMImerged(size(selected_PAJMImerged,2)) = toc(start);
end

delete(gcp('nocreate'));
