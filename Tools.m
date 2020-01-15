%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%merge the results
% path = '/home/hengl/matlab/bin/scripts/PAGreedy2.0/utils/';
% name = 'dbworld';
% 
% selected_PAUNImerged = [];
% time_PAUNImerged = [];
% window_trend = [];
% boot = [];
% 
% ind = 8;
% for i = 1:3
%     
%     path = strcat(name, num2str(i));
%     
%     temp = load(path, 'selected');
%     temp = temp.('selected');
%     temp = temp(1:ind,:);
%     selected_PAUNImerged = [selected_PAUNImerged; temp];
%     
%     temp = load(path, 'time');
%     temp = temp.('time');
%     temp = temp(1:ind,:);
%     time_PAUNImerged = [time_PAUNImerged; temp];
%     
%     temp = load(path, 'window');
%     temp = temp.('window');   
%     temp = temp(1:ind,:);
%     window_trend = [window_trend; temp];
%     
%     temp = load(path, 'bootstrap');
%     temp = temp.('bootstrap');
%     temp = temp(1:ind,:);
%     boot = [boot; temp];
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%measure the consistency
% temp1 = selected_UNIcached;
% temp2 = selected_PAUNImerged;
% 
% r = min([size(temp1,2) size(temp2,2)]);
% % Consistency
% index = [];
% for j=1:r
%     c = size(intersect(temp1(1:j),temp2(1:j)),2);
%     index = [index (c*f - j*j)/(j*f - j*j)];
% end
% mean(index)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The interfering zone displaying
fid = fopen('/home/hengl/matlab/bin/scripts/PAGreedy2.0/names2');
file_name = {};
for i=1:11
    line = fgets(fid);
    if ischar(line)
        file_name{i} = line;
    end
end
expression = '[A-Za-z0-9-_]*';
for i=1:11
    str = regexp(file_name{i}, expression);
    file_name{i} = file_name{i}(1: str(1,2)-2);
end


h = figure;
box on;
grid on;
algorithm = 'UNI';
for i = 1:11
    dataset = file_name{i};
    path = strcat('PAGreedy2.0/results/IZ/', algorithm, '/', dataset, '.mat');
    
    load(path, 'window_trend',  'f')
    selected_vector = load(path, strcat('selected_PA', algorithm, 'merged'));
    selected_vector = selected_vector.(strcat('selected_PA', algorithm, 'merged'));
    time_vector = load(path, strcat('time_PA', algorithm, 'merged'));
    time_vector = time_vector.(strcat('time_PA', algorithm, 'merged'));


%     load(path, 'window', 'selected', 'time', 'f')
%     time_vector = time;
%     selected_vector = selected;
%     window_trend = window;
    
    r = 0;
    ind = 0;
    for j = 1:size(selected_vector,1)
        temp = size(find(selected_vector(j,:) ~= 0),2);
        if r < temp
            r = temp;
            ind = j;
        end
    end
    time_vector = time_vector(ind,:);
    selected_trend = [];  
    r = 0;
    for j = time_vector
        r = r+1;
        if j ~= 0
            selected_trend = [selected_trend r];
        end
    end
    window_trend = window_trend(:, 1:size(selected_trend,2));
    if size(window_trend,1) > 1
        window_trend = mean(window_trend);
    end
    plot(selected_trend/f, window_trend/f, 'LineWidth', 2);
    hold on;
end
legend('dbworld', 'gisette', 'P53', 'drivFace', 'leukemia', 'duke-breast-cancer', 'arcene', 'dexter', 'TCGA-PANCAN', 'dorothea', 'PEMS-train', 'Location', 'best')
title(strcat('Interfering Zone with ', algorithm, ' objective'))
xlabel('Portion of features selected', 'FontSize', 20)
ylabel('Normalized Interfering Zone Size', 'FontSize', 16)
axis tight
%saveas(h, strcat(algorithm, '_IZ'), 'eps2c')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measure the time in individual figures, and here we only use one M option
% warmstart = [[0.05, 0.1, 0.04, 0.05, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01]
%     [0.01, 0.15, 0.15, 0.05, 0.02, 0.02, 0.12, 0.06, 0.03, 0.01, 0.01]
%     [0.08, 0.1, 0.1, 0.001, 0.02, 0.03, 0.02, 0.1, 0.01, 0.01, 0.005]];
% warmstart = warmstart';
% 
% h = figure;
% hold on
% 
% % specify the algorithm name
% track = 3;
% algorithms = {'JMI', 'MRMR', 'UNI'};
% algorithm = algorithms(track);
% 
% %k here is the index of the dataset
% data_index = 10; 
% names = {'dbworld', 'gisette', 'P53', 'drivFace', 'leukemia', 'duke_breast_cancer', 'arcene', 'dexter', 'TCGA_PANCAN', 'dorothea', 'PEMS_train'};
% file_name = names{data_index};
% 
% path1 = strcat('/home/hengl/matlab/bin/scripts/PAGreedy2.0/results/traditional/', algorithm, '/', file_name, '.mat') ;
% path2 = strcat('/home/hengl/matlab/bin/scripts/PAGreedy2.0/results/active/', algorithm, 'active/', file_name, '.mat');
% load(path1{1})
% load(path2{1})
% 
% path3 = strcat('time_', algorithm, 'cached');
% path4 = strcat('time_PA', algorithm, 'merged');
% temp1 = load(path1{1}, path3{1});
% temp1 = temp1.(path3{1});       %benchmark
% temp2 = load(path2{1}, path4{1});
% temp2 = temp2.(path4{1});       %PAGreedy
% 
% temp2(find(temp2 == 0)) = NaN;
% 
% r = min([size(temp1,2) size(temp2,2)]);
% 
% temp1 = temp1(1:r);
% temp2 = temp2(1:r);
% 
% witch = ceil(f*warmstart(data_index, track));
% 
% cheese = temp2;
% witch2 = witch;
% 
% if isnan(cheese(witch)) ~= 1
%     cheese(witch+1:size(cheese,2)) = cheese(witch+1:size(cheese,2)) - cheese(witch) + temp1(witch);
% else
%     for j = 1:witch
%         if (isnan(cheese(witch+j)) ~= 1) || (isnan(cheese(witch-j)) ~= 1)
%             if isnan(cheese(witch+j)) ~= 1
%                 witch2 = witch + j;
%             else
%                 witch2 = witch - j;
%             end
%             break;
%         end
%     end
%     cheese(witch2+1:size(cheese,2)) = cheese(witch2+1:size(cheese,2)) - cheese(witch2) + temp1(witch2);
% end
% 
% cheese(1:witch2) = temp1(1:witch2);
% temp2 = cheese;
% 
% 
% plot((1:r)/f, temp1/10000, 'LineWidth', 5);
% hold on
% 
% z = 1:r;
% grid on;
% box on;
% plot(z(find(isnan(temp2) ~= 1))/f, temp2(find(isnan(temp2) ~= 1))/10000, 'LineWidth', 5);
% 
% lgd = legend('traditional', 'PAGreedy', 'Location', 'best');
% set(lgd, 'FontSize', 25);
% xlabel('% of features selected', 'FontSize', 25)
% ylabel('Runtime(s) \times 10^4', 'FontSize', 25)
% axis tight
% % title(strcat(algorithm, '-', file_name), 'FontSize', 25);
% % saveas(h, strcat(algorithm, '_', file_name), 'eps2c');

% file = {'dbworld', 'gisette', 'p53', 'drivFace', 'leukemia', 'duke_breast_cancer', 'arcene', 'dexter', 'TCGA_PANCAN', 'dorothea', 'PEMS_train'};
% for i = 1:3
%     %line = file{i};
%     line = '';
%     for j = 1:11
%         line = strcat(line, '&', num2str(warmstart(i,j)));
%     end
%     line = strcat(line, '\\');
%     disp(line)
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the interfering zone for news20

% f = 1355191;
% x = 1:selected_trend(size(selected_trend, 2));
% y = zeros(1, selected_trend(size(selected_trend, 2)));
% for i = 1:size(selected_trend,2)
%     y(selected_trend(i)) = window_trend(i);
% end


%filling the gaps in the warmstarts for JMI
% y(1) = f;
% y(500) = 1354691;
% y(1000) = 1354191;
% y(1500) = 1353691;
% y(2000) = 1353191;
% y(2500) = 1352691;
% y(3000) = 301217;
% y(3500) = 118978;
% y(4000) = 78032;

% filling the gaps in the warmstarts for MRMR
% y(1) = f;
% y(500) = 1354553;
% y(1000) = 1353738;
% y(1500) = 1352653;
% y(2000) = 1350431;
% y(2500) = 1287747;
% y(3000) = 105130;
% y(3500) = 42333;

% filling the gaps in the warmstarts for CMI
% y(1000) = 33159;
% y(1500) = 19679;
% y(2000) = 22442;
% y(2500) = 28855;
% y(3000) = 30134;
% y(3500) = 45221;
% y(4000) = 78032;
% coeff = ones(1,2)/2;
% y = filter(coeff, 1, y);
% y(1) = f;
% y(500) = 1354524;


% h = figure;
% hold on
% plot(x(y ~= 0)/f, y(y~=0)/f, 'LineWidth', 2)
% 
% 
% xlabel('Portion of features selected', 'FontSize', 16)
% ylabel('Normalized Interfering Zone size', 'FontSize', 16)
% %title('Interfering Zone Size', 'FontSize', 15)
% axis tight
% box on
% grid on
% saveas(h, 'news20_CMI', 'eps2c');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot the time comparison for news20
% algorithm = 'UNI';
% file_name = 'news20';
% temp1 = time_UNI; %benchmark
% temp2 = time_PAUNImerged; %PAGreedy
% 
% r = min([size(temp1,2) size(temp2,2)]);
% 
% temp1 = temp1(1:r);
% temp2 = temp2(1:r);
% 
% h = figure;
% hold on
% plot((1:r)/f, temp1/10000, 'LineWidth', 5);
% 
% z = 1:r;
% grid on;
% box on;
% plot(z(find(temp2 ~= 0))/f, temp2(find(temp2 ~= 0))/10000, 'LineWidth', 5);
% 
% lgd = legend('traditional', 'PAGreedy', 'Location', 'best');
% set(lgd, 'FontSize', 25)
% xlabel('% of selected', 'FontSize', 25)
% ylabel('Runtime(s) \times 10^4', 'FontSize', 25)
% axis tight
%title(strcat(algorithm, '-', file_name), 'FontSize', 25);
% saveas(h, strcat(algorithm, '_', file_name), 'eps2c');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the time comaprison with multiple M's

% gisette = [[0.1, 0.001, 0.06];[0.1, 0.001, 0.06];[0.1, 0.001, 0.06];[0.1, 0.001, 0.06]];
% dexter = [[0.01, 0.07, 0.04]; [0.01, 0.07, 0.04]; [0.01, 0.07, 0.04]; [0.01, 0.07, 0.04]];
% pems = [[0.001, 0.005, 0.001]; [0.001, 0.005, 0.001]; [0.001, 0.005, 0.001]; [0.001, 0.005, 0.001]];
% 
% h = figure;
% hold on
% 
% % specify the algorithm name
% track = 2;
% algorithms = {'JMI', 'MRMR', 'UNI'};
% algorithm = algorithms(track);
% 
% %k here is the index of the dataset
% data_index = 3; 
% warmstart = pems;
% names = {'gisette', 'dexter', 'PEMS_train'};
% file_name = names{data_index};
% 
% path1 = strcat('/home/hengl/matlab/bin/scripts/PAGreedy2.0/results/traditional/', algorithm, '/', file_name, '.mat') ;
% path2 = strcat('/home/hengl/matlab/bin/scripts/PAGreedy2.0/results/active_VariSizeM/', algorithm, '/', file_name, '.mat');
% load(path1{1})
% load(path2{1})
% 
% path3 = strcat('time_', algorithm, 'cached');
% path4 = 'time';
% temp1 = load(path1{1}, path3{1});
% temp1 = temp1.(path3{1});       %benchmark
% temp2 = load(path2{1}, path4);
% temp2 = temp2.(path4);       %PAGreedy
% 
% 
% times = temp2;
% maximum = 0;
% 
% for m = 1:4
%     
%     temp2 = times(m ,:);
%     temp2(find(temp2 == 0)) = NaN;
% 
%     r = min([size(temp1,2) size(temp2,2)]);
% 
%     temp1 = temp1(1:r);
%     temp2 = temp2(1:r);
% 
%     witch = ceil(f*warmstart(m, track));
% 
%     cheese = temp2;
%     witch2 = witch;
% 
%     if isnan(cheese(witch)) ~= 1
%         cheese(witch+1:size(cheese,2)) = cheese(witch+1:size(cheese,2)) - cheese(witch) + temp1(witch);
%     else
%         for j = 1:witch
%             if (isnan(cheese(witch+j)) ~= 1) || (isnan(cheese(witch-j)) ~= 1)
%                 if isnan(cheese(witch+j)) ~= 1
%                     witch2 = witch + j;
%                 else
%                     witch2 = witch - j;
%                 end
%                 break;
%             end
%         end
%         cheese(witch2+1:size(cheese,2)) = cheese(witch2+1:size(cheese,2)) - cheese(witch2) + temp1(witch2);
%     end
% 
%     cheese(1:witch2) = temp1(1:witch2);
%     temp2 = cheese;
% 
%     z = 1:r;
%     grid on;
%     box on;
%     plot(z(find(isnan(temp2) ~= 1))/f, temp2(find(isnan(temp2) ~= 1))/10000, 'LineWidth', 3);
%     
%     if maximum < max(temp2)
%         maximum = max(temp2);
%     end
% end
% 
% %temp1 = temp1(1:1000);
% temp1 = temp1(find(temp1 <= 1.3*maximum));
% r = size(temp1, 2);
% 
% plot((1:r)/f, temp1(1:r)/10000, 'LineWidth', 3);
% lgd = legend('M=5', 'M=15', 'M=25', 'M=45', 'traditional', 'Location', 'best');
% set(lgd, 'FontSize', 15);
% xlabel('% of features selected', 'FontSize', 25)
% ylabel('Runtime(s) \times 10^4', 'FontSize', 25)
% axis tight
% 
