
trial = 1; % Change this to whatever simulation round you want to visualize


%% Segmented Data -- Comment This Section Out if you don't want this figure

%total_cells = ErrorPlusImage(trial).SegImage;
cell_points = zeros(size(total_cells));
cells = unique(total_cells);


for i = 1:length(cells)
    [row, col , z ] = ind2sub(size(total_cells),find(total_cells == i));
    
    s = 10 * ones(length(row),1); % size of edge pixels
    color = zeros(length(row),1);
    
    for j = 1:length(color)
        color(j) = i;
    end
    
    figure(10)
    scatter3(row(:),col(:),z(:),s(:,1),color(:,1),'filled')
    hold on

end

caxis([1,length(cells)]);
colormap jet;
title('Segmented Data')
xlabel('row')
ylabel('column')

%% Synthetic Data (non-convolved) -- Comment Out if you don't want this figure


% total_cells2 = ErrorPlusImage(trial).Labeled_SynthImage;
% cell_points = zeros(size(total_cells2));
% cells = unique(total_cells2);
% 
% 
% for i = 1:length(cells)
%     [row, col , z ] = ind2sub(size(total_cells2),find(total_cells2 == i));
%     
%     s = 10 * ones(length(row),1); % size of edge pixels
%     color = zeros(length(row),1);
%     
%     for j = 1:length(color)
%         color(j) = i;
%     end
%     
%     figure(11)
%     scatter3(row(:),col(:),z(:),s(:,1),color(:,1),'filled')
%     hold on
%     
% end
% 
% caxis([1,length(cells)]);
% colormap jet;
% title('Synthetic Cells')
% xlabel('row')
% ylabel('column')


%% Convolved Synthetic Data -- Comment OUt if you don't want this figure

% 
% total_cells3 = im_label2;
% cell_points = zeros(size(total_cells));
% cells = unique(total_cells);
% 
% 
% for i = 1:length(cells)
%     [row, col , z ] = ind2sub(size(total_cells3),find(total_cells3 == i));
%     
%     s = 10 * ones(length(row),1); % size of edge pixels
%     color = zeros(length(row),1);
%     
%     for j = 1:length(color)
%         color(j) = i;
%     end
%     
%     figure(12)
%     scatter3(row(:),col(:),z(:),s(:,1),color(:,1),'filled')
%     hold on
% 
% end
% 
% caxis([1,length(cells)]);
% colormap jet;
% title('Convolved Synthetic Data')
% xlabel('row')
% ylabel('column')