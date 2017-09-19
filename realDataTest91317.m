
field1 = 'Shape2D';
field2 = 'Concavity';
field3 = 'Cell';
TotalShape2 = struct(field1, [], field2, []);

for sim_round = [1:4,6:24]
    clearvars -except sim_round TotalShape2
    shape2D = [];
    clc
    field1 = 'ConcPoints';
    ConcStruct = struct(field1,[]);
    total_conc = [];
    
    
    dim  = 3;%input('Number of D''s (2/3) : ');
    ref_channel = 2; % Change to most in-focus channel. Probably 2/green or 3/blue
    ref_slice = 1;
    slices2D = 2; % How many frames above and below reference frame (e.g. 4 = reference frame +/- 4 frames)
    pix_size = .130; %Microns
    
    trials = 1; % Change this to change number of simulations run
    pixel_error = zeros(trials,4);
    int_thresh = .00001; % Intensity Threshold
    convolve_thresh = .05; % Threshold for Voxels to include in Convolved data
    ee_thresh = 1.5;  % <--- Splitting threshold, you can change this
    shape3D_thresh = 2;
    concavity_thresh = .3;
    background_thresh = .25;
    pix_neigh = 11; %floor((.08/pix_size)*12);
    volume_thresh = 500;
    slice_thresh = 3;
    zangle_thresh = 1;
    dist_thresh = 4;
    low_pass_check = 0; % 1 = on. Change to 0 if you want to turn it off
    gap_thresh = 5;
    
    % Path to main file (i.e. channel) that you will use for segmentation
    %     filepath = strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/August_24_17_convert/-SgrS(no_plasmid)/t=20/sample',num2str(cell_num)]);
    %     filepath_dic = strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/August_24_17_convert/-SgrS(no_plasmid)/t=20/dic',num2str(cell_num),'.tif']);
    %
    filepath = strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/cell', num2str(sim_round)]);
    
    
    [slice, stack_o, stack_red, stack_green, stack_blue, stack_back, slices, red_back, green_back, blue_back] = imFormat(filepath,ref_channel,dim,ref_slice,slices2D);
    
    %     for ig = 1:size(stack_o,3)
    %         stack_o(:,:,ig) = mat2gray(imread(filepath_dic));
    %     end
    
    % Uncomment these if you have the proper files described above
    se = [1 1 1; 1 1 1 ; 1 1 1]; % Structuring Element for basic Erosion
    % and dilation
    
    field1 = 'Stack_Number';
    field2 = 'Objects'; % All Objects, single and multi, labeled
    field3 = 'Center';
    field4 = 'Weighted_Center';
    field5 = 'Area';
    field6 = 'Ellipticity';
    field7 = 'Cell_Labels';
    field8 = 'Probability';
    field9 = 'Mask'; % Single Cell Selections
    field10 = 'All'; %All Objects, single and multi
    field11 = 'Original'; % Original Image
    field12 = 'Boundaries';
    field13 = 'Background';
    field14 = 'Non_single';
    part1 = struct(field1, [] , field2, [], field3, [], field4, [], field5, [], field6, [], field7, [], field8, [], field9, [], field10, [], field11, [], field12, [], field13, [], field14, []); % , field14, [], field15, []); %Table for Part 1
    
    %
    
    stack2 = zeros(size(stack_o));
    xdim = size(stack_o,1);
    ydim = size(stack_o,2);
    
    edge_cut = 2;
    
    for g = 1:slices
        %for g = 9
        strcat(['Working on Frame ' , num2str(g), ' ... '])
        stack2(:,:,g) = anisodiff2D(stack_o(:,:,g),3,1/7,30,1);
        I=stack_o(:,:,g);
        [r,c] = size(I);
        I2 = stack2(:,:,g);
        if dim == 2
            I_low_pass = low_pass(I2,.025);
        else
            I_low_pass = low_pass(I2,.05);
        end
        
        a(:,:,g) = imdilate(imerode(bradley(stack2(:,:,g),[pix_neigh,pix_neigh],int_thresh),se),se);
        if low_pass_check == 1
            I_LP = I2 - I_low_pass;
            a(:,:,g) = a(:,:,g) .* imdilate(im2bw(I_LP,.1),se);
        end
        b = bwlabel(a(:,:,g),4);
        a_temp = a(:,:,g);
        
        for i = 1:max(max(b))
            if  sum(sum(((b==i).*stack2(:,:,g))))/cellArea(b,i) < background_thresh
                a_temp(b==i) = 0;
                b(b == i) = 0;
            end
        end
        
        a(:,:,g) = a_temp;
        a(:,:,g) = imclearborder(a(:,:,g));
        a(1:edge_cut,:,g) = 0; a(r-edge_cut+1:r,:,g) = 0; a(:,1:edge_cut,g) = 0; a(:,c-edge_cut+1:c,g) = 0;
        BW = bwareaopen(a(:,:,g),10);
        [edges,BW] = edgeBreak(BW);
        BW = imfill(imclearborder(smallID(bwlabel(BW,4))),'holes');
        
        tic
        
        objects = bwlabel(BW,4);
        
        num = max(objects(:));
        
        ellipse_error = zeros(num,1);
        test_ellipse = {};
        
        area=zeros(num,1);
        
        for i=1:num
            
            area(i) = cellArea(objects,i,pix_size);
        end
        
        for i = 1:num
            
            [ellipse1,test1] = ellipseError(objects,i);
            if isempty(ellipse1) == 1 || isempty(test1) == 1
                ellipse_error(i) = ee_thresh+1;
                continue
            else
                test_ellipse(i) = test1;
                ellipse_error(i) = ellipseTest(ellipse1,test1,area(i),pix_size);
                
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
                %                     uni_obs = unique(objects2);
                %                     for numb = 1:length(unique(objects2))-1
                %                         objects2(objects2 == uni_obs(numb+1)) = numb;
                %                     end
            end
        end
        
        uni_obs = unique(objects);
        for numb = 1:length(unique(objects))-1
            objects(objects == uni_obs(numb+1)) = numb;
        end
        num = max(objects(:));
        
        clear centers area ellipticity
        
        centers=zeros(max(objects(:)),2);
        
        for i=1:num
            
            centers(i,:) = cellCenter(objects,i);
            
        end
        
        %Calculate Areas of connected objects. Not exactly just adding up
        %pixels. Also takes into account surrounding pixels
        area=zeros(num,1);
        
        for i=1:num
            
            area(i) = cellArea(objects,i,pix_size);
        end
        
        area = [area area];
        
        ellipticity = zeros(num,4);
        
        % Re-done Ellipticity Calculation
        for i=1:num
            %for i=5
            
            ellipticity(i,:) = cellEllipse(objects,i);
            
        end
        
        % Center of mass Calculation
        
        weighted_centers=zeros(num,2);
        
        for i=1:num
            weighted_centers(i,:) = weightedCenter(I,objects,i);
        end
        
        ellipticity = [ellipticity ellipticity];
        
        %probs = zeros(num,1);
        
        
        
        
        mask = BW;
        cell_labels = zeros(num,1);
        q = 1; %Non-Single Cells
        r = 1; % Single Cells
        %[edge_temp,bound_temp,con_peaks] = edgeOptimize(objects);
        %Single Cell Prediction
        for i = 1:max(max(objects))
            [edge_temp,bound_temp,con_peaks] = edgeOptimize(objects,i);
            %if probs(i) < 1E-4 || ellipse_error(i) >= ee_thresh %  Non - Single Cells hopefully
            if ellipse_error(i) >= ee_thresh | con_peaks>=3
                area(i,2) = NaN;
                ellipticity(i,5:8) = NaN;
                mask(objects == i) = 0;
                cell_labels(i) = 1000*(2*g)+q;
                q = q+1;
            else
                
                cell_labels(i) = 1000*(2*g-1) + r;
                r = r+1;
            end
        end
        
        
        %
        % Structured Array with our Data
        part1(g).Stack_Number = g;
        part1(g).Cell_Labels = cell_labels;
        part1(g).Area = area ;
        part1(g).Objects = objects;
        part1(g).Center = centers;
        part1(g).Weighted_Center = weighted_centers;
        part1(g).Ellipticity = ellipticity;
        part1(g).Probability = ellipse_error;
        part1(g).Mask = mask;
        part1(g).All = BW;
        part1(g).Original = I;
        part1(g).Boundaries = edges;
        part1(g).Background = I_low_pass;
        
        
        %     figure(g);imshow(part1(g).Object_ID)
        
        
    end
    
    for g = 1:slices
        objects = part1(g).Objects;
        num = max(objects(:));
        for i = 1:num
            shape2D = [shape2D part1(g).Probability(i)];
        end
        
        
        
        for i = 1:num
            [~,new_bounds,~] = edgeOptimize(objects,i);
             total_conc = [total_conc; new_bounds{1,1}(:,4)];
        end
    end
    
    TotalShape2(sim_round).Shape2D = shape2D;
    TotalShape2(sim_round).Concavity = total_conc;
    TotalShape2(sim_round).Cell = a;
    
end