%function [selected_PAUNImerged, time_PAUNImerged, window_trend] = PAUNILcachedmerged(dataindex, k, b, M, beta, gamma, selected_UNI, time_UNI, dist_feature_set, score_dist_features, reles)

%input: dataindex is the data, k is the amount of features selected, b is the bootstrapping indices, M is amount of worker available
%output: selected features and time vector

load('/home/hengl/PAGreedy/Experiments/PAGmerged/UNIcachednews201500.mat')
dataindex = 8;
k = 5001;
M = 17;
worker = 17;
beta = 0.5;
gamma = 0.5;

addpath /home/hengl/MIToolbox
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
hy = h(y);

%Initialization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_PAUNImerged = time_UNI;
selected_PAUNImerged = selected_UNI; %keep the selected feature subset
reles = reles';
rele = reles(:);

maximalrele = max(rele);
feature_set = 1:f; %To get rid of the all-zero columns


parpool(worker); %start the parallel pool and distributing the feature space
poolobj = gcp;
addAttachedFiles(poolobj,{'mi.m', 'cmi.m', 'h.m'});

previous_selected = []; %this keeps track of selected subset previous step

start = tic; %Start the clock

window_trend = [];
selected_trend = [];
%selecting k features%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while size(selected_PAUNImerged,2) <= k %select k or more features  
    %calculating the score of each candidate feature
    parfor i = 1:M  
        features = dist_feature_set(i,:);
        scores = score_dist_features(i,:); %old scores
        temp_scores = [];
        for p = features
            temp = 0;
            if p ~= 0 %we denote the selected/unavailable features as 0
                for q = previous_selected
                   temp = temp + beta*mi(full(data(:,p)), full(data(:,q))) - gamma*cmi(full(data(:,p)), full(data(:,q)), y);
                end
                temp_scores = [temp_scores temp];
            else
                temp_scores = [temp_scores 10000];
            end
        end
        score_dist_features(i,:) = scores - temp_scores;
    end
    
    %sort the candidate features and get the sequence
    scores = score_dist_features'; 
    scores = scores(:);
    [B, I] = sort(scores, 'descend');
    sorted_features = feature_set(I(B >= -10000));
    B = B(B >= -10000);

    %for each candidate calculate its boundary and verify the candidates in parallel 
    flag = ones(1, M+1); %1 means no intrusion
    intruders = zeros(1, M+1); %mark the intruders
    
    MIcorrespond = []; %prepare to store the temperory MI's 
    CMIcorrespond = []; %prepare to store the temperory CMI's 
    cell = {};          %the cached comparisons table
    
    windows = [];
    summation = [];
    
    parfor i = 2:M+1   %find the interfering zone
        boundary = 0;
        pre = []; %joint informativeness with features located ahead
        
        for s = 1:i-1 %find the boundary of the interfering zone 
            MI = mi(full(data(:,sorted_features(i))),full(data(:,sorted_features(s))));
            CMI = cmi(full(data(:,sorted_features(i))),full(data(:,sorted_features(s))),y);
            temp1 = beta * MI - gamma * CMI;
            temp2 = -1*gamma*rele(feature_set == sorted_features(s))+gamma*min(h(full(data(:, sorted_features(s)))), maximalrele);   
            pre = [pre temp1];
            boundary = boundary - (temp1 + temp2); 
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
        
        cell{i} = B(i) - sum(pre); 
    end
    window_trend = [window_trend max(windows)];
    disp(windows)
    windows = max(windows);
    

    parfor i = 1:M
        for j = i+1:windows
            temp = mi(full(data(:,sorted_features(i))),full(data(:,sorted_features(j))));
            MIcorrespond = [MIcorrespond; [sorted_features(i) sorted_features(j) temp]];
            temp = cmi(full(data(:,sorted_features(i))),full(data(:,sorted_features(j))),y);
            CMIcorrespond = [CMIcorrespond; [sorted_features(i) sorted_features(j) temp]];
        end
    end
    MImatrix = sparse([MIcorrespond(:,1);MIcorrespond(:,2)], [MIcorrespond(:,2);MIcorrespond(:,1)],[MIcorrespond(:,3);MIcorrespond(:,3)]); %store the MIs and make it symmetric
    CMImatrix = sparse([CMIcorrespond(:,1);CMIcorrespond(:,2)], [CMIcorrespond(:,2);CMIcorrespond(:,1)],[CMIcorrespond(:,3);CMIcorrespond(:,3)]); %store the useful CMIs and make it symmetric
    
    
    
    parfor i = 2:M+1   %identification
        temperory = [];
        intrude_max_objective = -10000; % to keep track of the score of the strongest intruder
        for p = i+1:windows        %iterate all possible intruders until the boundary reaches
            objective = B(i) - B(p) - summation(i);
            right = B(p);
            for q = 1:i-1  %each candidate is evaluated w.r.t. its predecessor 
                MI = MImatrix(sorted_features(p), sorted_features(q));
                CMI = CMImatrix(sorted_features(p), sorted_features(q));
                temp = beta*MI -gamma*CMI;
                right = right - temp;
                objective =  objective + temp;
            end
            temperory = [temperory; [cell{i} right]];
            if objective < 0 
                flag(i) = 0;
                if (B(i) - summation(i) - objective) > intrude_max_objective;
                    intruders(i) = sorted_features(p);
                    intrude_max_objective = B(i) - summation(i) - objective;
                end
            end 
        end 
        cell{i} = temperory;
    end
    disp(sorted_features(1:M+1))  %**testing**%
    disp(flag) %**testing**%
    disp(intruders)  %**testing**%
        
    
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
                        comparison(:,1) = comparison(:,1) + beta*MImatrix(originalfront, sorted_features(t)) - beta*MImatrix(originalrear, sorted_features(t));
                        comparison(:,1) = comparison(:,1) - gamma*CMImatrix(originalfront, sorted_features(t)) + gamma*CMImatrix(originalrear, sorted_features(t));
                        for z = 1:size(comparison,1)
                            if (z + t) == newplace
                                tempcell = cell{i};
                                comparison(z,2) = tempcell(1,1);
                                for w = i:t-1
                                    comparison(z,2) = comparison(z,2) - beta*MImatrix(sorted_features(w), originalfront) + gamma*CMImatrix(sorted_features(w), originalfront);
                                end
                            else
                                comparison(z,2) = comparison(z,2) + beta*MImatrix(originalfront, sorted_features(t+z)) - beta*MImatrix(originalrear, sorted_features(t+z));
                                comparison(z,2) = comparison(z,2) - gamma*CMImatrix(originalfront, sorted_features(t+z)) + gamma*CMImatrix(originalrear, sorted_features(t+z));
                            end
                        end
                    else
                        for m = i:newplace-1
                            for z = 1:size(comparison,1)
                                comparison(z,1) = comparison(z,1) - beta*MImatrix(originalfront, sorted_features(m)) + gamma*CMImatrix(originalfront, sorted_features(m));
                                comparison(z,2) = comparison(z,2) - beta*MImatrix(sorted_features(newplace+z), sorted_features(m)) + gamma*CMImatrix(sorted_features(newplace+z), sorted_features(m));
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
            scores(index) = -10001;
        end
        dist_feature_set(i,:) = features;
        score_dist_features(i,:) = scores;
    end
    selected_PAUNImerged = [selected_PAUNImerged previous_selected];
    selected_trend = [selected_trend size(selected_PAUNImerged,2)];
    
    disp('selecting...')
    disp(size(selected_PAUNImerged,2))

    
    time_PAUNImerged(size(selected_PAUNImerged,2)) = toc(start);
end

delete(gcp('nocreate'));

