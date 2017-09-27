clear all
clc


dim  = 3;%input('Number of D''s (2/3) : ');
ref_channel = 2; % Change to most in-focus channel. Probably 2/green or 3/blue
ref_slice = 5;
slices2D = 2; % How many frames above and below reference frame (e.g. 4 = reference frame +/- 4 frames)
pix_size = .130; %Microns

int_thresh = .00001; % Intensity Threshold
convolve_thresh = .05; % Threshold for Voxels to include in Convolved data
ee_thresh = 1.22;  % <--- Splitting threshold, you can change this
shape3D_thresh = 2;
concavity_thresh = .3;
background_thresh = .25;
pix_neigh = 11; %floor((.08/pix_size)*12); <- If you have no idea, try this
volume_thresh = 200;
slice_thresh = 3;
zangle_thresh = 1;
dist_thresh = 4;
low_pass_check = 0; % 1 = on. Change to 0 if you want to turn it off
gap_thresh = 5;

% Path to main file (i.e. channel) that you will use for segmentation
%filepath = strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/September_3_2017_convert/manX_gfp_no_plasmid/t20/sample',num2str(cell_num)]);
filepath = '/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/1_14_17_invasion_converted/0hr/wt_sopE_G1_crop2';


[slice, stack_o, stack_red, stack_green, stack_blue, stack_back, slices, red_back, green_back, blue_back] = imFormat(filepath,ref_channel,dim,ref_slice,slices2D);


se = [1 1 1; 1 1 1 ; 1 1 1]; % Structuring Element for basic Erosion and dilation

field1 = 'Stack_Number';
field2 = 'Objects'; % All Objects, single and multi, labeled
field3 = 'Cell_Labels';
field4 = 'Probability';
field5 = 'Mask'; % Single Cell Selections
field6 = 'All'; %All Objects, single and multi
field7 = 'Original'; % Original Image
field8 = 'Background';
field9 = 'Non_Single';
part1 = struct(field1, [] , field2, [], field3, [], field4, [], field5, [], field6, [], field7, [], field8, [], field9, []);  %Table for Part 1


stack2 = zeros(size(stack_o));
I3 = stack2;
xdim = size(stack_o,1);
ydim = size(stack_o,2);
edge_cut = 2;

for g = slice

    stack2(:,:,g) = anisodiff2D(stack_o(:,:,g),1,1/7,30,1);
    I=stack_o(:,:,g);
    I2 = stack2(:,:,g);
    if dim == 2
        I_low_pass = low_pass(I2,.025);
    else
        I_low_pass = low_pass(I2,.05);
    end
    I3(:,:,g) = stack2(:,:,g)-I_low_pass;
end

for g = slice
%for g = 18
    strcat(['Working on Frame ' , num2str(g), ' ... '])
    stack3 = I3(:,:,g)./max(I3(:));
    [r,c] = size(stack3);
    a(:,:,g) = imdilate(imerode(bradley(stack3,[pix_neigh,pix_neigh],int_thresh),se),se);
    b = bwlabel(a(:,:,g),4);
    a_temp = a(:,:,g);
    
    for i = 1:max(max(b))
        if  sum(sum(((b==i).*stack3)))/cellArea(b,i) < background_thresh
            a_temp(b==i) = 0;
            b(b == i) = 0;
        end
    end
    
    a(:,:,g) = a_temp;
    a(:,:,g) = imclearborder(a(:,:,g));
    a(1:edge_cut,:,g) = 0; a(r-edge_cut+1:r,:,g) = 0; a(:,1:edge_cut,g) = 0; a(:,c-edge_cut+1:c,g) = 0;
    BW = bwareaopen(a(:,:,g),10);
    [~,BW] = edgeBreak(BW);
    BW = imfill(imclearborder(smallID(bwlabel(BW,4))),'holes');
    
    objects = bwlabel(BW,4);
    num = max(objects(:));
    
    for i = 1:num
        [~,~,con_peaks] = edgeOptimize(objects,i);
        if con_peaks>=5
            objects(objects==i) = 0;
        end
    end
    
    ellipse_error = zeros(num,1);
    test_ellipse = {};
    
    for i = 1:num
        
        [ellipse1,test1] = ellipseError(objects,i);
        if isempty(ellipse1) == 1 || isempty(test1) == 1
            ellipse_error(i) = ee_thresh+1;
            continue
        else
            test_ellipse(i) = test1;
            ellipse_error(i) = ellipseTest(ellipse1,test1,cellArea(objects,i,pix_size),pix_size);
            
        end
    end
    
    for i = 1:num
        
        if ellipse_error(i) < ee_thresh
            objects(objects==i) = 0;
            object_temp = zeros(size(objects));
            for l = 1:length(test_ellipse{i})
                object_temp(test_ellipse{i}(l,1),test_ellipse{i}(l,2)) = i;
            end
            object_temp = imfill(object_temp);
            objects(object_temp == i) = i;
            objects = smallID(imfill(objects));
            
        end
    end
    
    uni_obs = unique(objects);
    for numb = 1:length(unique(objects))-1
        objects(objects == uni_obs(numb+1)) = numb;
    end
    num = max(objects(:));
    
    clear centers area ellipticity
    
    mask = BW;
    cell_labels = zeros(num,1);
    q = 1; %Non-Single Cells
    r = 1; % Single Cells
    non_single = [];
    
    %Single Cell Prediction
    for i = 1:num
        [~,~,con_peaks] = edgeOptimize(objects,i);

        if ellipse_error(i) >= ee_thresh || con_peaks>=3
            mask(objects == i) = 0;
            cell_labels(i) = 1000*(2*g)+q;
            q = q+1;
            non_single = [non_single i];
        else
            
            cell_labels(i) = 1000*(2*g-1) + r;
            r = r+1;
        end
    end
    
    
    %
    % Structured Array with our Data
    part1(g).Stack_Number = g;
    part1(g).Cell_Labels = cell_labels;
    part1(g).Objects = objects;
    part1(g).Probability = ellipse_error;
    part1(g).Mask = mask;
    part1(g).All = BW;
    part1(g).Original = stack_o(:,:,g);
    part1(g).Background = I_low_pass;
    part1(g).Non_Single = non_single;
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cluster Splitting

field1 = 'Stack_Number';
field2 = 'Objects'; % All Objects, single and multi, labeled
field3 = 'Center';
field4 = 'Area';
field5 = 'Ellipticity';
field6 = 'Cell_Labels';
field7 = 'Probability';
field8 = 'Mask'; % Single Cell Selections
field9 = 'All'; % All Objects, single and multi

part2 = struct(field1, [] , field2, [], field3, [], field4, [], field5, [], field6, [], field7, [], field8, [], field9, []); %Table for Part 2

split = cell(1,slices);
for g = slice
    %for g = 4
    strcat(['splitting frame ', num2str(g), ' ... '])
    
    objects = part1(g).Objects;
    objects2 = zeros(size(objects));
    non_single = part1(g).Non_Single;
    
    
    I2 = part1(g).All;
        
    for i = non_single
        
        if i > length(part1(g).Probability)
            continue
        end
        
        clear edge_temp bound_temp
        [split_im,~] = concave_split(objects,i,pix_size,ee_thresh);
        I2(objects == i) = 0;
        I2 = I2 + split_im;
    end
        
    objects2 = bwlabel(smallID(bwlabel(I2,4)),4);
    
    num = max(objects2(:));
    
    clear centers area ellipticity
    
    ellipse_error = zeros(num,1);
    test_ellipse = {};
   
    for round_id = 1:2
        for i = 1:num
            
            [ellipse1,test1] = ellipseError(objects2,i);
            
            if isempty(ellipse1) == 1 || isempty(test1) == 1
                ellipse_error(i) = ee_thresh+1;
                continue
            else
                test_ellipse(i) = test1;
                ellipse_error(i) = ellipseTest(ellipse1,test1,cellArea(objects2,i,pix_size),pix_size);
                
            end
        end
        
        for i = 1:num
            
            objects2(objects2==i) = 0;
            object_temp = zeros(size(objects2));
            for l = 1:length(test_ellipse{i})
                object_temp(test_ellipse{i}(l,1),test_ellipse{i}(l,2)) = i;
            end
            object_temp = imfill(object_temp);
            objects2(object_temp == i) = i;
            objects2 = smallID(imfill(objects2));
            
        end
        
        uni_obs = unique(objects2);
        for numb = 1:length(unique(objects2))-1
            objects2(objects2 == uni_obs(numb+1)) = numb;
        end
        num = max(objects2(:));
    end
    
    clear centers area ellipticity
    
  
    
    area=zeros(num,1);
    
    for i=1:num
        area(i) = cellArea(objects2,i,pix_size);
    end
    
    centers=zeros(max(objects2(:)),2);
    
    for i=1:num
        centers(i,:) = cellCenter(objects2,i);
    end
    
    ellipticity = zeros(num,4);
    
    % Re-done Ellipticity Calculation
    for i=1:num
        ellipticity(i,:) = cellEllipse(objects2,i);
    end
    
    ellipticity = [ellipticity ellipticity];
    
    
    cell_labels = zeros(num,1);
    q = 1; %Non-Single Cells
    r = 1; % Single Cells
    
    
    %Single Cell Prediction
    for i = 1:num
        if ellipse_error(i) >= ee_thresh %  Non - Single Cells hopefully
            ellipticity(i,5:8) = NaN;
            mask(objects2 == i) = 0;
            cell_labels(i) = 1000*(2*g)+q;
            q = q+1;
        else
            
            cell_labels(i) = 1000*(2*g-1) + r;
            r = r+1;
        end
    end
    
    part2(g).Stack_Number = g;
    part2(g).Cell_Labels = cell_labels;
    part2(g).Area = area ;
    part2(g).Objects = objects2;
    part2(g).Center = centers;
    part2(g).Ellipticity = ellipticity;
    part2(g).Probability = ellipse_error;
    part2(g).Mask = mask;
    part2(g).All = I2;
    clear cells
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Recombination

if dim == 3
    for g = 1:slices-1
        
        g
        
        if isempty(part2(g).Area) == 1 || isempty(part2(g+1).Area) == 1
            continue
        end
        
        for i= 1:length(part2(g).Area(:,1))    % i = Id'd cell (single or not) in Frame g
            
            if part2(g).Probability(i) > ee_thresh
                continue
            end
            
            clear ellipse4 ellipse6 possible_cell
            se_temp1 = strel('line',5,90-part2(g).Ellipticity(i,2));
            se_temp2 = strel('line',7,90-part2(g).Ellipticity(i,2));
            ob_temp = part2(g).Objects==i;
            [ellipse_temp,~] = ellipseError(part2(g).Objects,i);
            if isempty(ellipse_temp)
                continue
            end
            ob_dilate1 = imdilate(ob_temp,se_temp1);
            ob_erode1 = imerode(ob_temp,se_temp1);
            ob_dilate2 = imdilate(ob_temp,se_temp2);
            ob_erode2 = imerode(ob_temp,se_temp2);
            sum_erode1 = sum(sum(ob_erode1));
            sum_erode2 = sum(sum(ob_erode2));
            
            cell_distances = 1000*ones(length(part2(g+1).Area(:,1)),2);
            for j = 1:length(part2(g+1).Area(:,1))
                
                if part2(g+1).Probability(j) > ee_thresh
                    continue
                end
                
                cell_distances(j,1) = sqrt((part2(g).Center(i,1)-part2(g+1).Center(j,1))^2 + (part2(g).Center(i,2)-part2(g+1).Center(j,2))^2);
                cell_distances(j,2) = j;
                
            end
            [min_dist, possible_cell] = min(cell_distances(:,1));
            z_angle = atan(min_dist);
            if min_dist < dist_thresh
                
                [ellipse1,~] = ellipseError(part2(g).Objects,i);
                [ellipse2,~] = ellipseError(part2(g+1).Objects,possible_cell);
                [ellipse3,~] = ellipseError(ob_dilate1,1);
                if sum_erode1 >= 15
                    [ellipse4,~] = ellipseError(ob_erode1,1);
                else
                    ellipse4 = {};
                end
                [ellipse5,~] = ellipseError(ob_dilate2,1);
                if sum_erode2 >= 15
                    [ellipse6,~] = ellipseError(ob_erode2,1);
                else
                    ellipse6 = {};
                end
                
                deltaX = part2(g).Center(i,1) - part2(g+1).Center(possible_cell,1);
                deltaY = part2(g).Center(i,2) - part2(g+1).Center(possible_cell,2);
                ellipse1{1,1}(:,1) = ellipse1{1,1}(:,1) - round(deltaX);
                ellipse1{1,1}(:,2) = ellipse1{1,1}(:,2) - round(deltaY);
                ellipse3{1,1}(:,1) = ellipse3{1,1}(:,1) - round(deltaX);
                ellipse3{1,1}(:,2) = ellipse3{1,1}(:,2) - round(deltaY);
                if isempty(ellipse4) == 0
                    ellipse4{1,1}(:,1) = ellipse4{1,1}(:,1) - round(deltaX);
                    ellipse4{1,1}(:,2) = ellipse4{1,1}(:,2) - round(deltaY);
                end
                ellipse5{1,1}(:,1) = ellipse5{1,1}(:,1) - round(deltaX);
                ellipse5{1,1}(:,2) = ellipse5{1,1}(:,2) - round(deltaY);
                if isempty(ellipse6) == 0
                    ellipse6{1,1}(:,1) = ellipse6{1,1}(:,1) - round(deltaX);
                    ellipse6{1,1}(:,2) = ellipse6{1,1}(:,2) - round(deltaY);
                end
                
                if ellipseTest(ellipse1,ellipse2,100,pix_size) < shape3D_thresh || ellipseTest(ellipse3,ellipse2,100,pix_size) < shape3D_thresh || ellipseTest(ellipse4,ellipse2,100,pix_size) < shape3D_thresh || ellipseTest(ellipse5,ellipse2,100,pix_size) < shape3D_thresh || ellipseTest(ellipse6,ellipse2,100,pix_size) < shape3D_thresh   %|| (abs(part2(g).Ellipticity(i,2)-part2(g+1).Ellipticity(possible_cell,2))) < 7 + (sim_round/4)
                    possible_cell
                    part2(g+1).Cell_Labels(cell_distances(possible_cell,2)) = part2(g).Cell_Labels(i);
                elseif ellipseTest(ellipse1,ellipse2,100,pix_size) > shape3D_thresh && ellipseTest(ellipse3,ellipse2,100,pix_size) > shape3D_thresh && ellipseTest(ellipse4,ellipse2,100,pix_size) > shape3D_thresh && ellipseTest(ellipse5,ellipse2,100,pix_size) > shape3D_thresh && ellipseTest(ellipse6,ellipse2,100,pix_size) > shape3D_thresh && g <= slices - 2
                    'yes'
                    cell_distances2 = 1000*ones(length(part2(g+2).Area(:,1)),2);
                    for k = 1:length(part2(g+2).Area(:,1))
                        
                        if part2(g+2).Probability(k) > ee_thresh
                            continue
                        end
                        
                        cell_distances2(k,1) = sqrt((part2(g).Center(i,1)-part2(g+2).Center(k,1))^2 + (part2(g).Center(i,2)-part2(g+2).Center(k,2))^2);
                        cell_distances2(k,2) = k;
                        
                    end
                    [min_dist, possible_cell] = min(cell_distances2(:,1));
                    if min_dist < dist_thresh + 4
                        [ellipse1,test1] = ellipseError(part2(g).Objects,i);
                        [ellipse2,test2] = ellipseError(part2(g+2).Objects,possible_cell);
                        [ellipse3,~] = ellipseError(ob_dilate1,1);
                        if sum_erode1 >= 15
                            [ellipse4,~] = ellipseError(ob_erode1,1);
                        else
                            ellipse4 = {};
                        end
                        [ellipse5,~] = ellipseError(ob_dilate2,1);
                        if sum_erode2 >= 15
                            [ellipse6,~] = ellipseError(ob_erode2,1);
                        else
                            ellipse6 = {};
                        end
                        
                        deltaX = part2(g).Center(i,1) - part2(g+2).Center(possible_cell,1);
                        deltaY = part2(g).Center(i,2) - part2(g+2).Center(possible_cell,2);
                        ellipse1{1,1}(:,1) = ellipse1{1,1}(:,1) - round(deltaX);
                        ellipse1{1,1}(:,2) = ellipse1{1,1}(:,2) - round(deltaY);
                        ellipse3{1,1}(:,1) = ellipse3{1,1}(:,1) - round(deltaX);
                        ellipse3{1,1}(:,2) = ellipse3{1,1}(:,2) - round(deltaY);
                        if isempty(ellipse4) == 0
                            ellipse4{1,1}(:,1) = ellipse4{1,1}(:,1) - round(deltaX);
                            ellipse4{1,1}(:,2) = ellipse4{1,1}(:,2) - round(deltaY);
                        end
                        ellipse5{1,1}(:,1) = ellipse5{1,1}(:,1) - round(deltaX);
                        ellipse5{1,1}(:,2) = ellipse5{1,1}(:,2) - round(deltaY);
                        if isempty(ellipse6) == 0
                            ellipse6{1,1}(:,1) = ellipse6{1,1}(:,1) - round(deltaX);
                            ellipse6{1,1}(:,2) = ellipse6{1,1}(:,2) - round(deltaY);
                        end
                        if ellipseTest(ellipse1,ellipse2,100,pix_size) < shape3D_thresh || ellipseTest(ellipse3,ellipse2,100,pix_size) < shape3D_thresh || ellipseTest(ellipse4,ellipse2,100,pix_size) < shape3D_thresh || ellipseTest(ellipse5,ellipse2,100,pix_size) < shape3D_thresh || ellipseTest(ellipse6,ellipse2,100,pix_size) < shape3D_thresh   %|| (abs(part2(g).Ellipticity(i,2)-part2(g+1).Ellipticity(possible_cell,2))) < 7 + (sim_round/4)
                            i
                            possible_cell
                            part2(g+2).Cell_Labels(cell_distances2(possible_cell,2)) = part2(g).Cell_Labels(i);
                        elseif ellipseTest(ellipse1,ellipse2,100,pix_size) > shape3D_thresh && ellipseTest(ellipse3,ellipse2,100,pix_size) > shape3D_thresh && ellipseTest(ellipse4,ellipse2,100,pix_size) > shape3D_thresh && ellipseTest(ellipse5,ellipse2,100,pix_size) > shape3D_thresh && ellipseTest(ellipse6,ellipse2,100,pix_size) > shape3D_thresh && g <= slices - 3
                            'yes'
                            cell_distances2 = 1000*ones(length(part2(g+3).Area(:,1)),2);
                            for k = 1:length(part2(g+3).Area(:,1))
                                if part2(g+3).Probability(k) > ee_thresh
                                    continue
                                end
                                
                                cell_distances2(k,1) = sqrt((part2(g).Center(i,1)-part2(g+3).Center(k,1))^2 + (part2(g).Center(i,2)-part2(g+3).Center(k,2))^2);
                                cell_distances2(k,2) = k;
                                
                            end
                            [min_dist, possible_cell] = min(cell_distances2(:,1));
                            if min_dist < dist_thresh + 6
                                [ellipse1,test1] = ellipseError(part2(g).Objects,i);
                                [ellipse2,test2] = ellipseError(part2(g+3).Objects,possible_cell);
                                [ellipse3,~] = ellipseError(ob_dilate1,1);
                                if sum_erode1 >= 15
                                    [ellipse4,~] = ellipseError(ob_erode1,1);
                                else
                                    ellipse4 = {};
                                end
                                [ellipse5,~] = ellipseError(ob_dilate2,1);
                                if sum_erode2 >= 15
                                    [ellipse6,~] = ellipseError(ob_erode2,1);
                                else
                                    ellipse6 = {};
                                end
                                
                                deltaX = part2(g).Center(i,1) - part2(g+3).Center(possible_cell,1);
                                deltaY = part2(g).Center(i,2) - part2(g+3).Center(possible_cell,2);
                                ellipse1{1,1}(:,1) = ellipse1{1,1}(:,1) - round(deltaX);
                                ellipse1{1,1}(:,2) = ellipse1{1,1}(:,2) - round(deltaY);
                                ellipse3{1,1}(:,1) = ellipse3{1,1}(:,1) - round(deltaX);
                                ellipse3{1,1}(:,2) = ellipse3{1,1}(:,2) - round(deltaY);
                                if isempty(ellipse4) == 0
                                    ellipse4{1,1}(:,1) = ellipse4{1,1}(:,1) - round(deltaX);
                                    ellipse4{1,1}(:,2) = ellipse4{1,1}(:,2) - round(deltaY);
                                end
                                ellipse5{1,1}(:,1) = ellipse5{1,1}(:,1) - round(deltaX);
                                ellipse5{1,1}(:,2) = ellipse5{1,1}(:,2) - round(deltaY);
                                if isempty(ellipse6) == 0
                                    ellipse6{1,1}(:,1) = ellipse6{1,1}(:,1) - round(deltaX);
                                    ellipse6{1,1}(:,2) = ellipse6{1,1}(:,2) - round(deltaY);
                                end
                                if ellipseTest(ellipse1,ellipse2,100,pix_size) < shape3D_thresh || ellipseTest(ellipse3,ellipse2,100,pix_size) < shape3D_thresh || ellipseTest(ellipse4,ellipse2,100,pix_size) < shape3D_thresh || ellipseTest(ellipse5,ellipse2,100,pix_size) < shape3D_thresh || ellipseTest(ellipse6,ellipse2,100,pix_size) < shape3D_thresh   %|| (abs(part2(g).Ellipticity(i,2)-part2(g+1).Ellipticity(possible_cell,2))) < 7 + (sim_round/4)
                                    
                                    i
                                    possible_cell
                                    part2(g+3).Cell_Labels(cell_distances2(possible_cell,2)) = part2(g).Cell_Labels(i);
                                end
                            end
                        end
                    end
                end
            end
            
            
        end
    end
    
    total_cells = zeros(size(I,1),size(I,2),slices);
    
    for g = 1:slices
        total_cells(:,:,g) = part2(g).Objects;
        for xj = 1:max(max(part2(g).Objects))
            total_cells(total_cells ==xj) = part2(g).Cell_Labels(xj,1);
        end
    end
    
    indices = unique(total_cells(:));
    index = 0;
    
    
    for r = 2:length(indices)
        index = index + 1 ;
        total_cells(total_cells == indices(r)) = index;
        
    end
    
    for rk = 1:index
        slice_sum = 0;
        for g = 1:slices
            if sum(sum(total_cells(:,:,g)==rk)) > 0
                slice_sum = slice_sum + 1;
            end
        end
        if slice_sum <= 1
            total_cells(total_cells == rk) = 0;
        end
    end
    indices = unique(total_cells(:));
    index = 0;
    
    for r = 2:length(indices)
        index = index + 1 ;
        total_cells(total_cells == indices(r)) = index;
        
    end
    
    indices = 1:index;
    %% 3D Candidate Search
    
    field1 = 'cell_slices';
    field2 = 'delta_theta';
    field3 = 'delta_phi';
    field4 = 'center';
    field5 = 'dist_thresh';
    
    cell_angles = struct(field1,[],field2,[],field3, [],field4, [],field5, []);
    
    for gss =  1:length(indices)
        for g = 1:slices
            if sum(sum(total_cells(:,:,g)==gss)) > 0
                first = g;
                for gs = first:slices
                    if sum(sum(total_cells(:,:,gs)==gss)) ~= 0
                        last = gs;
                    end
                end
                
                cell_angles(gss).cell_slices = [first,last];
                break
            end
        end
    end
    
    for gss = 1:length(indices)
        if cell_angles(gss).cell_slices(1) >= slices
            continue
        end
        
        cell_angles(gss).center = zeros(cell_angles(gss).cell_slices(2)-cell_angles(gss).cell_slices(1)+1,2);
        theta = zeros(cell_angles(gss).cell_slices(2)-cell_angles(gss).cell_slices(1)+1,2);
        delta_theta = zeros(cell_angles(gss).cell_slices(2)-cell_angles(gss).cell_slices(1),1);
        delta_phi = zeros(cell_angles(gss).cell_slices(2)-cell_angles(gss).cell_slices(1),1);
        dist_thresh = zeros(cell_angles(gss).cell_slices(2)-cell_angles(gss).cell_slices(1),1);
        first = cell_angles(gss).cell_slices(1);
        last = cell_angles(gss).cell_slices(2);
        
        for g = cell_angles(gss).cell_slices(1):cell_angles(gss).cell_slices(2)
            cell_angles(gss).center(g-first+1,:) = cellCenter(total_cells(:,:,g),gss);
        end
        
        for gb = 1:length(theta)-1
            
            if sum(sum(total_cells(:,:,first+gb-1)==gss)) > 0 && sum(sum(total_cells(:,:,first+gb)==gss)) > 0
                theta(gb,1) = cell_angles(gss).center(gb+1,1)-cell_angles(gss).center(gb,1); %delta row
                theta(gb,2) = cell_angles(gss).center(gb+1,2)-cell_angles(gss).center(gb,2); %delta column
                dist_thresh(gb) = distance(cell_angles(gss).center(gb+1,:),cell_angles(gss).center(gb,:));
                delta_theta(gb) = atan(theta(gb,1)/theta(gb,2));
                delta_phi(gb) = atan(distance(cell_angles(gss).center(gb+1,:),cell_angles(gss).center(gb,:)));
            elseif cell_angles(gss).cell_slices(1) == cell_angles(gss).cell_slices(1)
                theta(gb,:) = [0,0];
                delta_theta(gb) = 0;
                delta_phi(gb) = 0;
            else
                theta(gb,:) = theta(gb-1,:);
                delta_theta(gb) = delta_theta(gb-1);
                delta_phi(gb) = delta_phi(gb-1);
            end
        end
        
        delta_theta(isnan(delta_theta)) = [];
        delta_phi(isnan(delta_phi)) = [];
        
        
        cell_angles(gss).delta_theta = mean(delta_theta);
        cell_angles(gss).delta_phi = mean(delta_phi);
        cell_angles(gss).dist_thresh = mean(dist_thresh);
    end
    
    %%
    candidates = {};
    candidate_index = 1;
    
    for i = 1:length(cell_angles)-1
        
        for j = i+1:length(cell_angles)
            
            if abs(cell_angles(i).delta_phi-cell_angles(j).delta_phi) < zangle_thresh
                
                first_last = cell_angles(i).center(length(cell_angles(i).cell_slices(2)),:); % Last slice of first cell
                last_first = cell_angles(j).center(1,:); % First slice of second cell
                slice_diff = cell_angles(j).cell_slices(1)-cell_angles(i).cell_slices(2);
                anglei = cellEllipse(total_cells(:,:,cell_angles(i).cell_slices(2)),i);
                anglej = cellEllipse(total_cells(:,:,cell_angles(j).cell_slices(1)),j);
                
                
                if slice_diff < gap_thresh && slice_diff > 0 && abs(anglei(2)-anglej(2)) < 50
                    delta_row_i = slice_diff*tan(cell_angles(i).delta_phi)*cos(cell_angles(i).delta_theta);
                    delta_col_i = slice_diff*tan(cell_angles(i).delta_phi)*sin(cell_angles(i).delta_theta);
                    delta_row_j = -slice_diff*tan(cell_angles(j).delta_phi)*cos(cell_angles(j).delta_theta);
                    delta_col_j = -slice_diff*tan(cell_angles(j).delta_phi)*sin(cell_angles(j).delta_theta);
                    
                    ix_coord = max([int8(round(first_last(1)+delta_row_i)),1]);
                    iy_coord = max([int8(round(first_last(2)+delta_col_i)),1]);
                    jx_coord = max([int8(round(last_first(1)+delta_row_j)),1]);
                    jy_coord = max([int8(round(last_first(2)+delta_col_j)),1]);
                    
                    ix_coord = min([ix_coord,xdim]);
                    iy_coord = min([iy_coord,ydim]);
                    jx_coord = min([jx_coord,xdim]);
                    jy_coord = min([jy_coord,ydim]);
                    
                    ix_neigh = [max([1,ix_coord-3]):min([ix_coord+3,xdim])];
                    iy_neigh = [max([1,iy_coord-3]):min([iy_coord+3,ydim])];
                    jx_neigh = [max([1,jx_coord-3]):min([jx_coord+3,xdim])];
                    jy_neigh = [max([1,jy_coord-3]):min([jy_coord+3,ydim])];
                    
                    i_overlap = 0;
                    j_overlap = 0;
                    
                    for k = 1:length(ix_neigh)
                        for l = 1:length(iy_neigh);
                            if total_cells(ix_neigh(k),iy_neigh(l),cell_angles(j).cell_slices(1)) == j
                                j_overlap = j_overlap + 1;
                            end
                        end
                    end
                    
                    for k = 1:length(jx_neigh)
                        for l = 1:length(jy_neigh)
                            if total_cells(jx_neigh(k),jy_neigh(l),cell_angles(i).cell_slices(2)) == i
                                i_overlap = i_overlap + 1;
                            end
                        end
                    end
                    
                    
                    if i_overlap >= 1 || j_overlap >= 1
                        candidates{candidate_index} = [i,j,0];
                        candidate_index = candidate_index + 1;
                    end
                    
                end
            end
        end
    end
    
    for i = 1:length(candidates)-1
        for j = i+1:length(candidates)
            if candidates{i}(1) == candidates{j}(1) %Two Cells trying to combine with bottom cell
                candidates{i}(3) = 1; % 1 means split bottom cell
                candidates{j}(3) = 1;
            elseif candidates{i}(2) == candidates{j}(2) % Two cells trying to combine with top cell
                candidates{i}(3) = 2; % 2 means split top cell
                candidates{j}(3) = 2;
            end
        end
    end
    
    for k = 1:length(candidates)
        i = candidates{k}(1);
        j = candidates{k}(2);
        if candidates{k}(3) == 0 % Combine as usual
            total_cells(total_cells == j) = i;
            if isnan(cell_angles(i).delta_theta)
                delta_row_i2 = 0;
                delta_col_i2 = 0;
            else
                delta_row_i2 = tan(cell_angles(i).delta_phi)*cos(cell_angles(i).delta_theta);
                delta_col_i2 = tan(cell_angles(i).delta_phi)*sin(cell_angles(i).delta_theta);
            end
            im_temp = total_cells(:,:,cell_angles(i).cell_slices(2)) == i;
            [xnow,ynow] = ind2sub(size(im_temp),find(im_temp));
            
            for gn = cell_angles(i).cell_slices(2)+1:cell_angles(j).cell_slices(1)-1
                xnow = xnow + delta_row_i2;
                xnow(xnow<1) = 1;
                ynow = ynow + delta_col_i2;
                ynow(ynow<1) = 1;
                for ind = 1:length(xnow)
                    total_cells(int8(round(xnow(ind))),int8(round(ynow(ind))),gn) = i;
                end
            end
        elseif candidates{k}(3) == 1 % Two Cells trying to combine with bottom cell
            
            if isnan(cell_angles(j).delta_theta)
                delta_row_i2 = 0;
                delta_col_i2 = 0;
            else
                delta_row_i2 = -tan(cell_angles(j).delta_phi)*cos(cell_angles(j).delta_theta);
                delta_col_i2 = -tan(cell_angles(j).delta_phi)*sin(cell_angles(j).delta_theta);
            end
            im_temp = total_cells(:,:,cell_angles(j).cell_slices(1)) == j;
            [xnow,ynow] = ind2sub(size(im_temp),find(im_temp));
            
            for gn = (cell_angles(j).cell_slices(1)-1):-1:cell_angles(i).cell_slices(1)
                xnow = xnow + delta_row_i2;
                xnow(xnow<1) = 1;
                ynow = ynow + delta_col_i2;
                ynow(ynow<1) = 1;
                for ind = 1:length(xnow)
                    total_cells(int8(round(xnow(ind))),int8(round(ynow(ind))),gn) = j;
                end
            end
        elseif candidates{k}(3) == 2 % Two cells trying to combine with top cell
            if isnan(cell_angles(i).delta_theta)
                delta_row_i2 = 0;
                delta_col_i2 = 0;
            else
                delta_row_i2 = tan(cell_angles(i).delta_phi)*cos(cell_angles(i).delta_theta);
                delta_col_i2 = tan(cell_angles(i).delta_phi)*sin(cell_angles(i).delta_theta);
            end
            im_temp = total_cells(:,:,cell_angles(j).cell_slices(2)) == i;
            [xnow,ynow] = ind2sub(size(im_temp),find(im_temp));
            
            for gn = cell_angles(i).cell_slices(1)+1:cell_angles(j).cell_slices(2)
                xnow = xnow + delta_row_i2;
                xnow(xnow<1) = 1;
                ynow = ynow + delta_col_i2;
                ynow(ynow<1) = 1;
                for ind = 1:length(xnow)
                    total_cells(int8(round(xnow(ind))),int8(round(ynow(ind))),gn) = i;
                end
            end
        end
    end
    
    
    for i = 1:length(cell_angles)
        for g = cell_angles(i).cell_slices(1):cell_angles(i).cell_slices(2)
            if sum(sum(total_cells(:,:,g)==i)) > 0
                last = g;
            elseif sum(sum(total_cells(:,:,g)==i)) == 0
                sli_diff = g - last;
                delta_row_i = sli_diff*tan(cell_angles(i).delta_phi)*cos(cell_angles(i).delta_theta);
                delta_col_i = sli_diff*tan(cell_angles(i).delta_phi)*sin(cell_angles(i).delta_theta);
                im_temp = total_cells(:,:,last) == i;
                [xnow,ynow] = ind2sub(size(im_temp),find(im_temp));
                xnow = xnow + delta_row_i;
                ynow = ynow + delta_col_i;
                xnow(xnow<1) = 1;
                ynow(ynow<1) = 1;
                
                for ind = 1:length(xnow)
                    total_cells(int8(round(xnow(ind))),int8(round(ynow(ind))),g) = i;
                end
            end
        end
    end
    
    
    % Min Volume Cutoff (vol_c < x). Set to 200 right now.
    for rk = 1:index
        
        slice_sum = 0;
        
        for g = 1:slices
            if sum(sum(total_cells(:,:,g)==rk)) > 0
                slice_sum = slice_sum + 1;
            end
        end
        
        [~,vol_c] = cellCenter(total_cells,rk);
        
        if length(vol_c) <= volume_thresh || slice_sum <= slice_thresh
            total_cells(total_cells == rk) = 0;
        end
        
        
    end
    
    
    indices = unique(total_cells(:));
    index = 0;
    
    for r = 2:length(indices)
        index = index + 1 ;
        total_cells(total_cells == indices(r)) = index;
    end
    
    indices = 1:index;
    
    total_cells_pre = total_cells;
    
    indices = unique(total_cells(:));
    index = 0;
    
    for r = 2:length(indices)
        index = index + 1 ;
        total_cells(total_cells == indices(r)) = index;
    end
    
    indices = 1:index;
    
    field1 = 'cell_slices';
    field2 = 'delta_theta';
    field3 = 'delta_phi';
    field4 = 'center';
    field5 = 'dist_thresh';
    
    cell_angles = struct(field1,[],field2,[],field3, [],field4, [],field5, []);
    
    for gss =  1:length(indices)
        for g = 1:slices
            if sum(sum(total_cells(:,:,g)==gss)) > 0
                first = g;
                for gs = first:slices
                    if sum(sum(total_cells(:,:,gs)==gss)) ~= 0
                        last = gs;
                    end
                end
                
                cell_angles(gss).cell_slices = [first,last];
                break
            end
        end
    end
    
    for gss = 1:length(indices)
        if cell_angles(gss).cell_slices(1) >= slices
            continue
        end
        
        cell_angles(gss).center = zeros(cell_angles(gss).cell_slices(2)-cell_angles(gss).cell_slices(1)+1,2);
        theta = zeros(cell_angles(gss).cell_slices(2)-cell_angles(gss).cell_slices(1)+1,2);
        delta_theta = zeros(cell_angles(gss).cell_slices(2)-cell_angles(gss).cell_slices(1),1);
        delta_phi = zeros(cell_angles(gss).cell_slices(2)-cell_angles(gss).cell_slices(1),1);
        dist_thresh = zeros(cell_angles(gss).cell_slices(2)-cell_angles(gss).cell_slices(1),1);
        first = cell_angles(gss).cell_slices(1);
        last = cell_angles(gss).cell_slices(2);
        
        if gss == 5
            gss;
        end
        
        
        for g = cell_angles(gss).cell_slices(1):cell_angles(gss).cell_slices(2)
            cell_angles(gss).center(g-first+1,:) = cellCenter(total_cells(:,:,g),gss);
            
        end
        
        for gb = 1:length(theta)-1
            
            if sum(sum(total_cells(:,:,first+gb-1)==gss)) > 0 && sum(sum(total_cells(:,:,first+gb)==gss)) > 0
                theta(gb,1) = cell_angles(gss).center(gb+1,1)-cell_angles(gss).center(gb,1); %delta row
                theta(gb,2) = cell_angles(gss).center(gb+1,2)-cell_angles(gss).center(gb,2); %delta column
                dist_thresh(gb) = distance(cell_angles(gss).center(gb+1,:),cell_angles(gss).center(gb,:));
                delta_theta(gb) = atan(theta(gb,1)/theta(gb,2));
                delta_phi(gb) = atan(distance(cell_angles(gss).center(gb+1,:),cell_angles(gss).center(gb,:)));
            elseif cell_angles(gss).cell_slices(1) == cell_angles(gss).cell_slices(1)
                theta(gb,:) = [0,0];
                delta_theta(gb) = 0;
                delta_phi(gb) = 0;
            else
                theta(gb,:) = theta(gb-1,:);
                delta_theta(gb) = delta_theta(gb-1);
                delta_phi(gb) = delta_phi(gb-1);
            end
        end
        
        delta_theta(isnan(delta_theta)) = [];
        delta_phi(isnan(delta_phi)) = [];
        
        
        cell_angles(gss).delta_theta = mean(delta_theta);
        cell_angles(gss).delta_phi = delta_phi;
        cell_angles(gss).dist_thresh = mean(dist_thresh);
    end
    
elseif dim == 2
    total_cells = zeros(size(I,1),size(I,2),(2*slices2D+1));
    sliced = (slice-slices2D):(slice+slices2D);
    
    for g = 1:length(sliced)
        total_cells(:,:,g) = part2(slice).Objects;
    end
    
    indices = unique(total_cells(:));
    index = 0;
    
    for r = 2:length(indices)
        index = index + 1 ;
        total_cells(total_cells == indices(r)) = index;
    end
end


%% Final Characterizaton (Volume, Integrated Intensity, etc.)

field1 = 'Cell';
field2 = 'Volume';
field3 = 'Center';
field4 = 'Intensity';
field5 = 'Intensity_Red';
field6 = 'Intensity_Green';
field7 = 'Intensity_Blue';
field8 = 'Background';
field9 = 'Red_Background';
field9 = 'Green_Background';
field10 = 'Blue_Background';

part4 = struct(field1,[],field2,[],field3,[],field4,[],field5,[],field6,[],field7,[],field8,[],field9,[],field10,[]);


for i = 1:index
    i
    part4(i).Cell = i;
    [center, volume] = cellCenter(total_cells,i);
    part4(i).Volume = length(volume);
    part4(i).Center = center;
    part4(i).Intensity = cellIntensity(total_cells,stack2,i);
    part4(i).Background = cellIntensity(total_cells,stack_back,i);
    
    if exist('stack_red') == 1
        part4(i).Intensity_Red = cellIntensity(total_cells,stack_red,i);
        part4(i).Red_Background = cellIntensity(total_cells,red_back,i);
    end
    
    if exist('stack_blue') == 1
        part4(i).Intensity_Blue = cellIntensity(total_cells,stack_blue,i);
        part4(i).Blue_Background = cellIntensity(total_cells,blue_back,i);
    end
    
    if exist('stack_green') == 1
        part4(i).Intensity_Green = cellIntensity(total_cells,stack_green,i);
        part4(i).Green_Background = cellIntensity(total_cells,green_back,i);
    end
    
end


