%function [selected_PAMRMRmerged, time_PAMRMRmerged, window_trend] = PAMRMRLcachedmerged(dataindex, k, b, M, selected_MRMR, time_MRMR, dist_feature_set, score_dist_features, reles)

%This algorithm is also taking: 'selected_MRMR', 'time_MRMR', 'dist_feature_set', 'score_dist_features', 'reles' as inputs
%input: dataindex is the data, k is the amount of features selected, b is the bootstrapping indices, M is amount of worker available
%output: selected features and time vector

addpath /home/hengl/MIToolbox
load('/home/hengl/PAGreedy/Experiments/Greedy/MRMRcachednews201500.mat')
dataindex = 8;
k = 2501;
M = 48;
worker = 48;

%Data preparation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('/home/hengl/PAGreedy2.0/names2');
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
path = strcat('/home/hengl/PAGreedy2.0/data/', file_name{dataindex}, '.mat');
file = load(path);
data = struct2cell(file);
data = data{1};
f = size(data,2)-1;
y = full(data(:, size(data,2)));


%Initialization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_PAMRMRmerged = time_MRMR;
selected_PAMRMRmerged = selected_MRMR; %keep the selected feature subset

rele = reles(1:f);
feature_set = 1:f; %To get rid of the all-zero columns

sum_entropy = 0;  %this is the sum of all selected features' entropy, used in the bounds
for current = selected_PAMRMRmerged
    sum_entropy = sum_entropy + h(full(data(:, current))); 
end

parpool(worker); %start the parallel pool and distributing the feature space
poolobj = gcp;
addAttachedFiles(poolobj,{'mi.m', 'joint.m'});

previous_selected = []; %this keeps track of selected subset previous step

start = tic; %Start the clock

window_trend = [];
selected_trend = [];
%selecting k features%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while size(selected_PAMRMRmerged,2) <= k %select k or more features  
    %calculating the score of each candidate feature
    cursize = size(selected_PAMRMRmerged,2); %current selected subset size
    presize = size(previous_selected,2); %previous selected subset size
    if isempty(previous_selected) ~= 1
	    parfor i = 1:M  
		features = dist_feature_set(i,:);
		scores = score_dist_features(i,:); %old scores
		temp_scores = [];
		for p = features
		    temp = 0;
		    if p ~= 0 %we denote the selected/unavailable features as 0
		        temp = rele(feature_set == p) + (scores(features == p) - rele(feature_set == p))*(cursize - presize)/cursize;
		        for q = previous_selected
		           temp = temp - mi(full(data(:,p)),full(data(:,q)))/cursize;
		        end
		        temp_scores = [temp_scores temp];
		    else
		        temp_scores = [temp_scores -10000];
		    end
		end
		score_dist_features(i,:) = temp_scores;
	    end
    end
    %sort the candidate features and get the sequence
    scores = score_dist_features'; 
    scores = scores(:);
    [B, I] = sort(scores, 'descend');
    sorted_features = feature_set(I(B > -10000));
    B = B(B > -10000);

    %for each candidate calculate its boundary and verify the candidates in parallel 
    flag = ones(1, M+1); %1 means no intrusion
    intruders = zeros(1, M+1); %mark the intruders
    
    MIcorrespond = [];        %prepare to store the temperory MI's 
    cell = {};                %the cached comparisons
    windows = [];
    summation = [];
    
    parfor i = 2:M+1   %find the interfering zone
        boundary = 0;
        pre = []; %joint informativeness with features located ahead
        
        for s = 1:i-1 %find the boundary of the interfering zone 
            temp = mi(full(data(:,sorted_features(i))),full(data(:,sorted_features(s))));
            pre = [pre temp];
            boundary = boundary + temp/(cursize + i - 1);  
        end
        
        summation(i) = sum(pre)/(cursize + i - 1);
        windowsize = 1;
        
        %mapping the score to new space
        Bi_prime = rele(feature_set == sorted_features(i)) + (B(i) - rele(feature_set == sorted_features(i)))*cursize/(cursize + i - 1);
        for p = i+1:size(sorted_features,2)
            if Bi_prime - B(p) - boundary - sum_entropy/cursize/(cursize + i - 1) < 0
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
        
        cell{i} = Bi_prime - sum(pre)/(cursize + i - 1); 
    end
    
    disp(windows)
    windows = min((max(windows)+M), size(sorted_features,2));
    if windows > 25000
        windows = floor(mean(window_trend))
    end
    window_trend = [window_trend windows];
    
    
    parfor i = 1:M
        for j = i+1:windows
            temp = mi(full(data(:,sorted_features(i))),full(data(:,sorted_features(j))));
            MIcorrespond = [MIcorrespond; [sorted_features(i) sorted_features(j) temp]];
        end
    end
    MImatrix = sparse([MIcorrespond(:,1);MIcorrespond(:,2)], [MIcorrespond(:,2);MIcorrespond(:,1)],[MIcorrespond(:,3);MIcorrespond(:,3)]); %store the useful MI's and make it symmetric     
    
    
    parfor i = 2:M+1   %identification
        temperory = [];
        intrude_max_objective = -10000; % to keep track of the score of the strongest intruder
        for p = i+1:windows        %iterate all possible intruders until the boundary reaches
            Bi_prime = rele(feature_set == sorted_features(i)) + (B(i) - rele(feature_set == sorted_features(i)))*cursize/(cursize + i - 1);
            Bp_prime = rele(feature_set == sorted_features(p)) + (B(p) - rele(feature_set == sorted_features(p)))*cursize/(cursize + i - 1);
            objective = Bi_prime - Bp_prime - summation(i);
            right = Bp_prime;
            for q = 1:i-1  %each candidate is evaluated w.r.t. its predecessor 
                temp = MImatrix(sorted_features(p), sorted_features(q));
                right = right - temp/(cursize+i-1);
                objective =  objective + temp/(cursize+i-1);
            end
            temperory = [temperory; [cell{i} right]];
            if objective < 0 
                flag(i) = 0;
                if Bi_prime - summation(i) - objective > intrude_max_objective;
                    intruders(i) = sorted_features(p);
                    intrude_max_objective = Bi_prime - summation(i) - objective;
                end
            end 
        end 
        cell{i} = temperory;
    end
    %disp(sorted_features(1:M+1))  %**testing**%
    disp(flag) %**testing**%
    %disp(intruders)  %**testing**%      
    
    
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
                flag(i) = 1; flag(i+1:newplace) = 0; intruders(i+1:newplace) = 0; %Initialize identifications and Maximal Intruders
                previous_selected = [previous_selected sorted_features(i)];
                %disp(sorted_features(1:M+1))  %**testing**%
                for t = i+1:newplace  %updating comparisons
                    comparison = cell{t};
                    if t < newplace
                        comparison(:,1) = comparison(:,1) + (MImatrix(originalfront, sorted_features(t)) - MImatrix(originalrear, sorted_features(t)))/(cursize + t - 1);
                        for z = 1:size(comparison,1)
                            if (z + t) == newplace
                                tempcell = cell{i};
                                temper = tempcell(1,1);
                                comparison(z,2) = rele(feature_set == sorted_features(newplace)) + (temper - rele(feature_set == sorted_features(newplace)))*(cursize + i - 1)/(cursize + t - 1);
                                for w = i:t-1
                                    comparison(z,2) = comparison(z,2) - MImatrix(sorted_features(w), originalfront)/(cursize + t - 1);
                                end
                            else
                                comparison(z,2) = comparison(z,2) + (MImatrix(originalfront, sorted_features(t+z)) - MImatrix(originalrear, sorted_features(t+z)))/(cursize + t - 1);
                            end
                        end
                    else
                        comparison(:,1) = rele(feature_set == sorted_features(newplace)) + (comparison(:,1) - rele(feature_set == sorted_features(newplace)))*(cursize + i - 1)/(cursize + t - 1);
                        for m = i:newplace-1
                            comparison(:,1) = comparison(:,1) - MImatrix(originalfront, sorted_features(m))/(cursize + t - 1);
                            for z = 1:size(comparison,1)
                                if m == i %in case multiple mapping
                                    comparison(z,2) = rele(feature_set == sorted_features(newplace+z)) + (comparison(z,2) - rele(feature_set == sorted_features(newplace+z)))*(cursize + i - 1)/(cursize + t - 1);
                                end
                                comparison(z,2) = comparison(z,2) - MImatrix(sorted_features(newplace+z), sorted_features(m))/(cursize + t - 1);
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
                %disp(intruders)  %**testing**%
            else
                break;
            end
        end       
    end   
    index = find(flag == 0);
    for i = 1:index
        sum_entropy = sum_entropy + h(full(data(:, sorted_features(i))));
    end
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
            scores(index) = -10001;
        end
        dist_feature_set(i,:) = features;
        score_dist_features(i,:) = scores;
    end
    selected_PAMRMRmerged = [selected_PAMRMRmerged previous_selected];
    selected_trend = [selected_trend size(selected_PAMRMRmerged,2)];
    
    disp('selecting...')
    disp(size(selected_PAMRMRmerged,2))

    
    time_PAMRMRmerged(size(selected_PAMRMRmerged,2)) = toc(start) + time_MRMR(size(time_MRMR, 2));
end

delete(gcp('nocreate'));

