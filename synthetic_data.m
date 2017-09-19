%Single_volume = zeros(1,50);
load psf.mat
int_thresh = .1; % Intensity Threshold
convolve_thresh = .05; % Threshold for Voxels to include in Convolved data
ee_thresh = 1.5;  % <--- Splitting threshold, you can change this
for trial = 1;

    
    %close all
    
    % global  p4 p5 p6 a4 a5 a6 e3 e4 ConSet;
    % global TRb;
    n = 1; % # Cells
    field1 = 'cell';
    field2 = 'coordinates';
    synth = struct(field1,[],field2, []);
    slices = 20;
    xdim = 100;
    ydim = 100;
    
    im_s = zeros(xdim,ydim,slices);
    im_label = zeros(xdim,ydim,slices);
    
    
    for m = 1:n
        
        p4=int32(xdim*rand); p5=int32(ydim*rand); p6=int32(slices*rand);          % position (x,y,z) Surf 2 Green
        
        a4= rand + 2; a5 = a4; a6= 6 + 14*rand;           % superellipsoid 2 semi-axes Green
        
        e3=4; e4=2;                   % Cell Shape : e3 = 4 , e4 = 2 (Roundness)
        
        
        phi = 2*pi*rand;
        theta = pi*rand;
        psi = 2*pi*rand;
        
        euler1B = cos(.5*(phi-psi))*sin(.5*theta);
        euler2B = sin(.5*(phi-psi))*sin(.5*theta);
        euler3B = sin(.5*(phi+psi))*cos(.5*theta);   % Euler param. Surf 2 (orientation)
        
        ConSet='normals_dot';
        'surfaces'          % surfaces only       a)
        'normals_cross'     % n1 x n2             b)
        'normals_cross+n1'  % n1 x n2 + n2.r12    c)
        'normals_cross+n2'  % n1 x n2 + n2.r12    d)
        'normals_dot'       % n1.n2               e)
        'normals_dot+n1'    % n1.n2   + n1.r12    f)
        'normals_dot+n2'    % n1.n2   + n2.r12    g)
        'dist_n1_n2'        % n1.r12  + n2.r12    h)
        
       
        euler0B=sqrt(1.0 - euler1B^2 - euler2B^2 - euler3B^2);
       
        e3(1)=2.0/e3(1);  e4(1)=2.0/e4(1);
        
        
        
        TRb(1,1)= euler0B^2 + euler1B^2 - euler2B^2 - euler3B^2; % 1st column = TRA(1)
        TRb(2,1)= 2*(euler1B*euler2B + euler0B*euler3B);         %            = TRA(2)
        TRb(3,1)= 2*(euler1B*euler3B - euler0B*euler2B);         %            = TRA(3)
        TRb(1,2)= 2*(euler1B*euler2B - euler0B*euler3B);         % 2nd column = TRA(4)
        TRb(2,2)= euler0B^2 - euler1B^2 + euler2B^2 - euler3B^2; %            = TRA(5)
        TRb(3,2)= 2*(euler2B*euler3B + euler0B*euler1B);         %            = TRA(6)
        TRb(1,3)= 2*(euler1B*euler3B + euler0B*euler2B);         % 3rd column = TRA(7)
        TRb(2,3)= 2*(euler2B*euler3B - euler0B*euler1B);         %            = TRA(8)
        TRb(3,3)= euler0B^2 - euler1B^2 - euler2B^2 + euler3B^2; %            = TRA(9)
        
        figure1 = figure('PaperSize',[20.98 29.68],'Color',[1 1 1],...
            'NumberTitle','off','ToolBar','figure','MenuBar','none','Units','pixels');
        
        % Create axes
        axes1 = axes('Parent',gcf,'DataAspectRatio',[1 1 1],'TickDir','in',...
            'FontName','Calibri','FontWeight','light','FontSize',12);
        
        set(gcf,'color','white', 'Position',[20 50 1000 750],...
            'nextplot','replacechildren');
        view([120 20]);
        hold('all');
        lighting GOURAUD % FLAT, GOURAUD, PHONG, NONE
        camlight headlight
        
        nxy = 200; nxz = nxy;
        
        
        % Surface 2
        i=1;
        for w=-pi : 2*pi/nxy : pi;
            j=1;
            for n=-pi/2 : pi/nxz : pi/2;
                xB(j,i)=a4*sign(cos(n))*abs(cos(n)).^e3*sign(cos(w))*abs(cos(w)).^e4;
                yB(j,i)=a5*sign(cos(n))*abs(cos(n)).^e3*sign(sin(w))*abs(sin(w)).^e4;
                zB(j,i)=a6*sign(sin(n))*abs(sin(n)).^e3;
                j=j+1;
            end
            i=i+1;
        end
        
        
        xB_global= double(p4) + xB*TRb(1,1) + yB*TRb(1,2) + zB*TRb(1,3);
        yB_global= double(p5) + xB*TRb(2,1) + yB*TRb(2,2) + zB*TRb(2,3);
        zB_global= double(p6) + xB*TRb(3,1) + yB*TRb(3,2) + zB*TRb(3,3);
        
        x = int32(xB_global);
        y = int32(yB_global);
        z = int32 (zB_global);
        
        x = x - min(x(:)) + 1;
        y = y - min(y(:)) + 1;
        z = z - min(z(:)) + 1;
        
        surfaceB = surf(x,y,z,'Parent',gca,...
            'EdgeLighting','gouraud',...
            'FaceLighting','gouraud',...
            'LineWidth',1.5,...
            'FaceColor',[0.6 1.0 0.6],...
            'FaceAlpha',0.6,...
            'EdgeColor',[0.2 0.2 0.2]);
        
        axis equal; axis on; grid on;
        
        
        se = zeros(3,3,3);
        
        se(:,:,1) = [1,1,1;1,1,1;1,1,1];
        se(:,:,2) = [1,1,1;1,1,1;1,1,1];
        se(:,:,3) = [1,1,1;1,1,1;1,1,1];
        
        cellmat = zeros(max(x(:)),max(y(:)),max(z(:)));
        
        coords = [x(:),y(:),z(:)];
        coords = unique(coords,'rows');
        
        for i = 1:length(coords)
            cellmat(sub2ind(size(cellmat),coords(i,1),coords(i,2),coords(i,3))) = 1;
        end
        
        cellmat = imerode(imdilate(cellmat,se),se);
        
        for i = 1:size(cellmat,3)
            cellmat(:,:,i) = imfill(cellmat(:,:,i),'holes');
        end
        
        cellmat = m*bwlabeln(cellmat);
        
        [CellX,CellY,CellZ] = size(cellmat);
        
        
        cell_x = [p4-floor(size(cellmat,1)/2),p4+floor(size(cellmat,1)/2)];
        cell_y = [p5-floor(size(cellmat,2)/2),p5+floor(size(cellmat,2)/2)];
        cell_z = [p6-floor(size(cellmat,3)/2),p6+floor(size(cellmat,3)/2)];
        
        cell1x = cell_x(1) - 1;
        cell2x = xdim-cell_x(2);
        cell1y = cell_y(1)-1;
        cell2y = ydim-cell_y(2);
        cell1z = cell_z(1)-1;
        cell2z = slices-cell_z(2);
        
        if cell1x < 0
            cell_x(1) = 1;
            cellmat(1:abs(cell1x),:,:) = [];
            
            [CellX,CellY,CellZ] = size(cellmat);
        end
        
        if cell2x < 0
            cell_x(2) = xdim;
            cellmat(((CellX+cell2x):CellX),:,:) = [];
            
            [CellX,CellY,CellZ] = size(cellmat);
        end
        
        if cell1y < 0
            cell_y(1) = 1;
            cellmat(:,1:abs(cell1y),:) = [];
            
            [CellX,CellY,CellZ] = size(cellmat);
        end
        
        if cell2y < 0
            cell_y(2) = ydim;
            cellmat(:,(CellY+cell2y):CellY,:) = [];
            
            [CellX,CellY,CellZ] = size(cellmat);
        end
        
        if cell1z < 0
            cell_z(1) = 1;
            cellmat(:,:,1:abs(cell1z)) = [];
            [CellX,CellY,CellZ] = size(cellmat);
        end
        
        if cell2z < 0
            cell_z(2) = slices;
            cellmat(:,:,(CellZ+cell2z):CellZ) = [];
            
            [CellX,CellY,CellZ] = size(cellmat);
        end
        
        
        im_label_m = zeros(xdim,ydim,slices);
        im_label_m(cell_x(1):cell_x(1)+CellX-1,cell_y(1):cell_y(1)+CellY-1,cell_z(1):cell_z(1)+CellZ-1) = cellmat;
        
        
        C = convn(cellmat,psf);
        C = (.5*rand+.5)*(C./max(C(:)));
        
        [Cx,Cy,Cz] = size(C);
        
        
        coordinates_x = [p4-floor(size(C,1)/2),p4+floor(size(C,1)/2)];
        coordinates_y = [p5-floor(size(C,2)/2),p5+floor(size(C,2)/2)];
        coordinates_z = [p6-floor(size(C,3)/2),p6+floor(size(C,3)/2)];
        
        while sum(im_label_m.*im_label)>0
            if p4 > xdim/2
                cell_x = cell_x - 1;
                coordinates_x(1) = coordinates_x(1)-1;
                coordinates_x(2) = coordinates_x(2)-1;
            else
                cell_x = cell_x + 1;
                coordinates_x(1) = coordinates_x(1)+1;
                coordinates_x(2) = coordinates_x(2)+1;
            end
            
            if p5 > ydim/2
                cell_y = cell_y - 1;
                coordinates_y(1) = coordinates_y(1)-1;
                coordinates_y(2) = coordinates_y(2)-1;
            else
                cell_y = cell_y + 1;
                coordinates_y(1) = coordinates_y(1)+1;
                coordinates_y(2) = coordinates_y(2)+1;
            end
            
            if p6 > zdim/2
                cell_z = cell_z - 1;
                coordinates_z(1) = coordinates_z(1)-1;
                coordinates_z(2) = coordinates_z(2)-1;
            else
                cell_z = cell_z + 1;
                coordinates_z(1) = coordinates_z(1)+1;
                coordinates_z(2) = coordinates_z(2)+1;
            end
            
            im_label_m = zeros(xdim,ydim,slices);
            im_label_m(cell_x(1):cell_x(1)+CellX-1,cell_y(1):cell_y(1)+CellY-1,cell_z(1):cell_z(1)+CellZ-1) = cellmat;
            
            
        end
        
        c1x = coordinates_x(1) - 1;
        c2x = xdim-coordinates_x(2);
        c1y = coordinates_y(1)-1;
        c2y = ydim-coordinates_y(2);
        c1z = coordinates_z(1)-1;
        c2z = slices-coordinates_z(2);
        
        if c1x < 0
            coordinates_x(1) = 1;
            C(1:abs(c1x),:,:) = [];
            [Cx,Cy,Cz] = size(C);
        end
        
        if c2x < 0
            coordinates_x(2) = xdim;
            C((Cx+c2x):Cx,:,:) = [];
            [Cx,Cy,Cz] = size(C);
        end
        
        if c1y < 0
            coordinates_y(1) = 1;
            C(:,1:abs(c1y),:) = [];
            [Cx,Cy,Cz] = size(C);
        end
        
        if c2y < 0
            coordinates_y(2) = ydim;
            C(:,(Cy+c2y):Cy,:) = [];
            [Cx,Cy,Cz] = size(C);
        end
        
        if c1z < 0
            coordinates_z(1) = 1;
            C(:,:,1:abs(c1z)) = [];
            [Cx,Cy,Cz] = size(C);
            
        end
        
        if c2z < 0
            coordinates_z(2) = slices;
            C(:,:,(Cz+c2z):Cz) = [];
            [Cx,Cy,Cz] = size(C);
        end
        
        coordinates = [coordinates_x,coordinates_y,coordinates_z];
        
        k = 1;
        
        im_label(cell_x(1):cell_x(1)+CellX-1,cell_y(1):cell_y(1)+CellY-1,cell_z(1):cell_z(1)+CellZ-1) = cellmat;
        im_s(coordinates(1):coordinates(1)+size(C,1)-1,coordinates(3):coordinates(3)+size(C,2)-1,coordinates(5):coordinates(5)+size(C,3)-1) = im_s(coordinates(1):coordinates(1)+size(C,1)-1,coordinates(3):coordinates(3)+size(C,2)-1,coordinates(5):coordinates(5)+size(C,3)-1) + C;
        
    end
    
    im_s = im_s./max(im_s(:));
    
    close all
    se = [1 1 1; 1 1 1 ; 1 1 1]; % Structuring Element for basic Erosion
    % and dilation
    
    pix_size = .130; % in microns, e.g. if pixel size = 130 nm --> pix_size =.130
    pix_neigh = floor((.08/pix_size)*12);
    
    %%%%%%%%%%%%%%%%% -----------
    
    
    
    dim = 3;
    
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
    part1 = struct(field1, [] , field2, [], field3, [], field4, [], field5, [], field6, [], field7, [], field8, [], field9, [], field10, [], field11, [], field12, [], field13, []); % , field14, [], field15, []); %Table for Part 1
    
    %
    
    stack2 = anisodiff3D(im_s,1,3/44,30,1,[1,1,9]);
    
    
    [maskLabel,mask,ncells] = segmentation_threshold_Jing(stack2,400);
    lower_thresh = 1; upper_thresh = 100000; % Size Thresholds
    [ncells,maskLabel3,maskLabel3props,mask3] = segmentation_cellfilter(maskLabel,lower_thresh,upper_thresh,ncells);
    
    
    edge_cut = 2;
    
    for g = 1:slices
        %for g = 9
        strcat(['Working on Frame ' , num2str(g), ' ... '])
        I=im_s(:,:,g);
        [r,c] = size(I);
        I2 = stack2(:,:,g);
%         if dim == 2
%             I_low_pass = low_pass(I2,.025);
%         else
%             I_low_pass = low_pass(I2,.05);
%         end
%         stack_back(:,:,g) = I_low_pass;
%         I_LP = I2 - I_low_pass;
        a(:,:,g) = imfill(imdilate(imerode(bradley(stack2(:,:,g),[pix_neigh,pix_neigh],int_thresh),se),se),'holes');
%         a(:,:,g) = a(:,:,g) .* imdilate(im2bw(I_LP,.1),se);
        a(:,:,g) = imclearborder(a(:,:,g));
        a(1:edge_cut,:,g) = 0; a(r-edge_cut+1:r,:,g) = 0; a(:,1:edge_cut,g) = 0; a(:,c-edge_cut+1:c,g) = 0;
        BW = bwareaopen(a(:,:,g),5);
        [edges,BW] = edgeBreak(BW);
        BW = imclearborder(smallID(BW));
        
        tic
        
        
        objects = bwlabel(BW,4);
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
        
        % Probability Calculation. Compare to Single-cell PDF, which should
        % probably be refined at this point. Talk to Jingyi
        %     for i = 1:num
        %         probs(i) = cell_prob(ellipticity(i,1),area(i,1),ellip_compare,area_compare,density);
        %     end
        
        ellipse_error = zeros(num,1);
        
        for i = 1:num
            
            [ellipse1,test1] = ellipseError(objects,i);
            if isempty(ellipse1) == 1 || isempty(test1) == 1
                ellipse_error(i) = 1;
            else
                ellipse_error(i) = ellipseTest(ellipse1,test1,area(i),pix_size);
            end
        end
        
        mask = BW;
        cell_labels = zeros(num,1);
        q = 1; %Non-Single Cells
        r = 1; % Single Cells
        
        
        %Single Cell Prediction
        for i = 1:max(max(objects))
            %if probs(i) < 1E-4 || ellipse_error(i) >= ee_thresh %  Non - Single Cells hopefully
            if ellipse_error(i) >= ee_thresh
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
        
        
        
        %     figure(g);imshow(part1(g).Object_ID)
        
        
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Cluster Splitting
    
    field1 = 'Stack_Number';
    field2 = 'Objects'; % All Objects, single and multi, labeled
    field3 = 'Center';
    field4 = 'Weighted_Center';
    field5 = 'Area';
    field6 = 'Ellipticity';
    field7 = 'Cell_Labels';
    field8 = 'Probability';
    field9 = 'Mask'; % Single Cell Selections
    field10 = 'All'; % All Objects, single and multi
    field11 = 'Boundaries';
    part2 = struct(field1, [] , field2, [], field3, [], field4, [], field5, [], field6, [], field7, [], field8, [], field9, [], field10, [], field11, []); % , field14, [], field15, []); %Table for Part 1
    
    data_compare = zeros(xdim,ydim,slices);
    
    
    split = cell(1,slices);
    for g = 1:slices
        %for g = 7
        strcat(['splitting frame ', num2str(g), ' ... '])
        
        objects = ones(size(part1(g).All));
        objects2 = zeros(size(objects));
        split_round = 0;
        splits = {};
        while sum(sum(objects ~= objects2)) > 0  && split_round < 1
            
            split_round = split_round + 1;
            
            if split_round == 1
                objects = part1(g).Objects;
                cells = length(part1(g).Cell_Labels);
                I2 = part1(g).All;
                splits{1,g} = cell(1,cells);
            else
                objects = part2(g).Objects;
                cells = length(part2(g).Cell_Labels);
                I2 = part2(g).All;
            end
            
            
            
            for i = 1:cells
                if  split_round == 1
                    if i > length(part1(g).Probability)
                        continue
                    end
                    
                    clear edge_temp bound_temp
                    
                    [edge_temp,bound_temp,con_peaks] = edgeOptimize(objects);
                    
                    if (part1(g).Probability(i) > ee_thresh && max(bound_temp{i,1}(:,4)) > .2) || con_peaks(i,1) > 0 %Multi-Cell
                        i
                        [split_im,split_lines] = concave_split(objects,i,pix_size,ee_thresh);
                        
                        for new_i = 1:max(max(bwlabel(split_im,4)))
                            
                            split_im_temp = bwlabel(split_im,4);
                            
                            [ellipse1,test1] = ellipseError(split_im_temp,new_i);
                            if isempty(ellipse1) == 1 || isempty(test1) == 1
                                ellipse_error(new_i) = 2;
                            else
                                ellipse_error(new_i) = ellipseTest(ellipse1,test1,100,pix_size);
                                
                                [edge_temp,bound_temp,con_peaks] = edgeOptimize(split_im_temp);
                                
                                if (ellipse_error(new_i) > ee_thresh && max(bound_temp{new_i,1}(:,4)) > .2) || con_peaks(new_i,1) > 0
                                    [split_im1,~] = concave_split(split_im_temp,new_i,pix_size,ee_thresh);
                                    split_im(split_im_temp == new_i) = 0;
                                    split_im = split_im + split_im1;
                                end
                            end
                            
                            
                        end
                        
                        I2(objects == i) = 0;
                        I2 = I2 + split_im;
                        split{1,g}{1,i} = split_lines;
                    elseif mod(floor(part1(g).Cell_Labels(i)/100),2) == 0
                        i;
                        split{1,g}{1,i} = [];
                    end
                    
                elseif split_round > 1
                    if i > length(part1(g).Probability)
                        continue
                    end
                    
                    clear edge_temp bound_temp
                    
                    [edge_temp,bound_temp,con_peaks] = edgeOptimize(objects);
                    
                    if (part1(g).Probability(i) > ee_thresh && max(bound_temp{i,1}(:,4)) > .2) || con_peaks(i,1) > 0
                        [split_im,~] = concave_split(objects,i,pix_size,ee_thresh);
                        
                        for new_i = 1:max(max(bwlabel(split_im,4)))
                            
                            split_im_temp = bwlabel(split_im,4);
                            
                            [ellipse1,test1] = ellipseError(split_im_temp,new_i);
                            if isempty(ellipse1) == 1 || isempty(test1) == 1
                                ellipse_error(new_i) = 2;
                            else
                                ellipse_error(new_i) = ellipseTest(ellipse1,test1,100,pix_size);
                                
                                clear edge_temp bound_temp
                                
                                [edge_temp,bound_temp,con_peaks] = edgeOptimize(split_im_temp);
                                
                                if (ellipse_error(new_i) > ee_thresh && max(bound_temp{new_i,1}(:,4)) > .2) || con_peaks(new_i,1) > 0
                                    [split_im1,~] = concave_split(split_im_temp,new_i,pix_size,ee_thresh);
                                    split_im(split_im_temp == new_i) = 0;
                                    split_im = split_im + split_im1;
                                end
                            end
                            
                            
                        end
                        
                        I2(objects == i) = 0;
                        I2 = I2 + split_im;
                        %splits = {splits,split_lines};
                    end
                end
            end
            
            objects2 = bwlabel(smallID(I2),4);
            
            num = max(objects2(:));
            
            clear centers area ellipticity
            
            centers=zeros(max(objects2(:)),2);
            
            for i=1:num
                
                centers(i,:) = cellCenter(objects2,i);
                
            end
            
            %Calculate Areas of connected objects. Not exactly just adding up
            %pixels. Also takes into account surrounding pixels
            area=zeros(num,1);
            
            for i=1:num
                
                area(i) = cellArea(objects2,i,pix_size);
            end
            
            area = [area area];
            
            ellipticity = zeros(num,4);
            
            % Re-done Ellipticity Calculation
            for i=1:num
                %for i=5
                
                ellipticity(i,:) = cellEllipse(objects2,i);
                
            end
            
            % Center of mass Calculation
            
            weighted_centers=zeros(num,2);
            
            for i=1:num
                weighted_centers(i,:) = weightedCenter(part1(g).Original,objects2,i);
            end
            
            ellipticity = [ellipticity ellipticity];
            
            %probs = zeros(num,1);
            
            % Probability Calculation. Compare to Single-cell PDF, which should
            % probably be refined at this point. Talk to Jingyi
            %         for i = 1:num
            %             probs(i) = cell_prob(ellipticity(i,1),area(i,1),ellip_compare,area_compare,density);
            %         end
            
            ellipse_error = zeros(num,1);
            
            for i = 1:num
                
                [ellipse1,test1] = ellipseError(objects2,i);
                if isempty(ellipse1) == 1 || isempty(test1) == 1
                    ellipse_error(i) = 1;
                else
                    ellipse_error(i) = ellipseTest(ellipse1,test1,area(i),pix_size);
                end
                
                
            end
            
            mask = smallID(I2);
            cell_labels = zeros(num,1);
            q = 1; %Non-Single Cells
            r = 1; % Single Cells
            
            
            %Single Cell Prediction
            for i = 1:max(max(objects2))
                if ellipse_error(i) >= ee_thresh %  Non - Single Cells hopefully
                    area(i,2) = NaN;
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
            part2(g).Weighted_Center = weighted_centers;
            part2(g).Ellipticity = ellipticity;
            part2(g).Probability = ellipse_error;
            part2(g).Mask = mask;
            part2(g).All = I2;
            part2(g).Boundaries = edges;
            
        end
        
        
        clear cells
        
        data_compare(:,:,g) = I2;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Recombination
    
    
    for g = 1:slices-1
        %for g = 6
        g
        for i= 1:length(part2(g).Area(:,1))    % i = Id'd cell (single or not) in Frame g
            %for i = 2
            clear ellipse4 ellipse6
            se_temp1 = strel('line',5,90-part2(g).Ellipticity(i,2));
            se_temp2 = strel('line',7,90-part2(g).Ellipticity(i,2));
            ob_temp = part2(g).Objects==i;
            ob_dilate1 = imdilate(ob_temp,se_temp1);
            ob_erode1 = imerode(ob_temp,se_temp1);
            ob_dilate2 = imdilate(ob_temp,se_temp2);
            ob_erode2 = imerode(ob_temp,se_temp2);
            sum_erode1 = sum(sum(ob_erode1));
            sum_erode2 = sum(sum(ob_erode2));
            
            cell_distances = 1000*ones(length(part2(g+1).Area(:,1)),2);
            for j = 1:length(part2(g+1).Area(:,1))
                %if mod(floor(part3(g+1).Cell_Labels(j)/100),2) == 1  % single cell in frame g+1
                
                cell_distances(j,1) = sqrt((part2(g).Center(i,1)-part2(g+1).Center(j,1))^2 + (part2(g).Center(i,2)-part2(g+1).Center(j,2))^2);
                cell_distances(j,2) = j;
                %end
            end
            [min_dist, possible_cell] = min(cell_distances(:,1));
            z_angle = atan(min_dist);
            if min_dist < 5 %& (abs(part2(g).Ellipticity(i,2)-part2(g+1).Ellipticity(possible_cell,2)) < 10 || (part2(g).Ellipticity(i,1) < .82 && part2(g+1).Ellipticity(possible_cell,1) < .82))
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
                
                if ellipseTest(ellipse1,ellipse2,100,pix_size) < ee_thresh || ellipseTest(ellipse3,ellipse2,100,pix_size) < ee_thresh || ellipseTest(ellipse4,ellipse2,100,pix_size) < ee_thresh || ellipseTest(ellipse5,ellipse2,100,pix_size) < ee_thresh || ellipseTest(ellipse6,ellipse2,100,pix_size) < ee_thresh   %|| (abs(part2(g).Ellipticity(i,2)-part2(g+1).Ellipticity(possible_cell,2))) < 7 + (sim_round/4)
                    possible_cell
                    part2(g+1).Cell_Labels(cell_distances(possible_cell,2)) = part2(g).Cell_Labels(i);
                elseif ellipseTest(ellipse1,ellipse2,100,pix_size) > ee_thresh && ellipseTest(ellipse3,ellipse2,100,pix_size) > ee_thresh && ellipseTest(ellipse4,ellipse2,100,pix_size) > ee_thresh && ellipseTest(ellipse5,ellipse2,100,pix_size) > ee_thresh && ellipseTest(ellipse6,ellipse2,100,pix_size) > ee_thresh && g <= slices - 2
                    'yes'
                    cell_distances2 = 1000*ones(length(part2(g+2).Area(:,1)),2);
                    for k = 1:length(part2(g+2).Area(:,1))
                        %if mod(floor(part3(g+1).Cell_Labels(j)/100),2) == 1  % single cell in frame g+1
                        cell_distances2(k,1) = sqrt((part2(g).Center(i,1)-part2(g+2).Center(k,1))^2 + (part2(g).Center(i,2)-part2(g+2).Center(k,2))^2);
                        cell_distances2(k,2) = k;
                        %end
                    end
                    [min_dist, possible_cell] = min(cell_distances2(:,1));
                    if min_dist < 8 %& (abs(part2(g).Ellipticity(i,2)-part2(g+1).Ellipticity(possible_cell,2)) < 10 || (part2(g).Ellipticity(i,1) < .82 && part2(g+1).Ellipticity(possible_cell,1) < .82))
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
                        if ellipseTest(ellipse1,ellipse2,100,pix_size) < ee_thresh || ellipseTest(ellipse3,ellipse2,100,pix_size) < ee_thresh || ellipseTest(ellipse4,ellipse2,100,pix_size) < ee_thresh || ellipseTest(ellipse5,ellipse2,100,pix_size) < ee_thresh || ellipseTest(ellipse6,ellipse2,100,pix_size) < ee_thresh   %|| (abs(part2(g).Ellipticity(i,2)-part2(g+1).Ellipticity(possible_cell,2))) < 7 + (sim_round/4)
                            i
                            possible_cell
                            part2(g+2).Cell_Labels(cell_distances2(possible_cell,2)) = part2(g).Cell_Labels(i);
                        elseif ellipseTest(ellipse1,ellipse2,100,pix_size) > ee_thresh && ellipseTest(ellipse3,ellipse2,100,pix_size) > ee_thresh && ellipseTest(ellipse4,ellipse2,100,pix_size) > ee_thresh && ellipseTest(ellipse5,ellipse2,100,pix_size) > ee_thresh && ellipseTest(ellipse6,ellipse2,100,pix_size) > ee_thresh && g <= slices - 3
                            'yes'
                            cell_distances2 = 1000*ones(length(part2(g+3).Area(:,1)),2);
                            for k = 1:length(part2(g+3).Area(:,1))
                                %if mod(floor(part3(g+1).Cell_Labels(j)/100),2) == 1  % single cell in frame g+1
                                cell_distances2(k,1) = sqrt((part2(g).Center(i,1)-part2(g+3).Center(k,1))^2 + (part2(g).Center(i,2)-part2(g+3).Center(k,2))^2);
                                cell_distances2(k,2) = k;
                                %end
                            end
                            [min_dist, possible_cell] = min(cell_distances2(:,1));
                            if min_dist < 10 %& (abs(part2(g).Ellipticity(i,2)-part2(g+1).Ellipticity(possible_cell,2)) < 10 || (part2(g).Ellipticity(i,1) < .82 && part2(g+1).Ellipticity(possible_cell,1) < .82))
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
                                if ellipseTest(ellipse1,ellipse2,100,pix_size) < ee_thresh || ellipseTest(ellipse3,ellipse2,100,pix_size) < ee_thresh || ellipseTest(ellipse4,ellipse2,100,pix_size) < ee_thresh || ellipseTest(ellipse5,ellipse2,100,pix_size) < ee_thresh || ellipseTest(ellipse6,ellipse2,100,pix_size) < ee_thresh   %|| (abs(part2(g).Ellipticity(i,2)-part2(g+1).Ellipticity(possible_cell,2))) < 7 + (sim_round/4)
                                    
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
    
    
    
    
    % for rk = 1:index
    %     slice_sum = 0;
    %     cell_slice = [];
    %     for g = 1:slices
    %         if sum(sum(total_cells(:,:,g)==rk)) > 0
    %             slice_sum = slice_sum + 1;
    %             cell_slice = [cell_slice, g];
    %         end
    %
    %     end
    %
    %     if slice_sum <= 3
    %         total_cells(total_cells == rk) = 0;
    %     end
    %
    % end
    
    
    indices = unique(total_cells(:));
    index = 0;
    
    for r = 2:length(indices)
        index = index + 1 ;
        total_cells(total_cells == indices(r)) = index;
    end
    
    indices = 1:index;
    %% 3D Candidate Search
    
%     field1 = 'cell_slices';
%     field2 = 'delta_theta';
%     field3 = 'delta_phi';
%     field4 = 'center';
%     
%     cell_angles = struct(field1,[],field2,[],field3, [],field4, []);
%     
%     for gss =  1:length(indices)
%         for g = 1:slices
%             if sum(sum(total_cells(:,:,g)==gss)) > 0
%                 first = g;
%                 for gs = first:slices
%                     if sum(sum(total_cells(:,:,gs)==gss)) ~= 0
%                         last = gs;
%                     end
%                 end
%                 
%                 cell_angles(gss).cell_slices = [first,last];
%                 break
%             end
%         end
%     end
%     
%     for gss = 1:length(indices)
%         cell_angles(gss).center = zeros(cell_angles(gss).cell_slices(2)-cell_angles(gss).cell_slices(1)+1,2);
%         theta = zeros(cell_angles(gss).cell_slices(2)-cell_angles(gss).cell_slices(1)+1,2);
%         delta_theta = zeros(cell_angles(gss).cell_slices(2)-cell_angles(gss).cell_slices(1),1);
%         delta_phi = zeros(cell_angles(gss).cell_slices(2)-cell_angles(gss).cell_slices(1),1);
%         first = cell_angles(gss).cell_slices(1);
%         last = cell_angles(gss).cell_slices(2);
%         
%         if gss == 5
%             gss;
%         end
%         
%         
%         for g = cell_angles(gss).cell_slices(1):cell_angles(gss).cell_slices(2)
%             cell_angles(gss).center(g-first+1,:) = cellCenter(total_cells(:,:,g),gss);
%             %         theta(g-first+1,1) = cell_angles(gss).center(g-first+1,1)/cell_angles(gss).center(g-first+1,2));
%             %         theta(g-first+1,1) = atan(cell_angles(gss).center(g-first+1,1)/cell_angles(gss).center(g-first+1,2));
%             %
%         end
%         
%         for gb = 1:length(theta)-1
%             if sum(sum(total_cells(:,:,first+gb-1)==gss)) > 0 && sum(sum(total_cells(:,:,first+gb)==gss)) > 0
%                 theta(gb,1) = cell_angles(gss).center(gb+1,1)-cell_angles(gss).center(gb,1); %delta row
%                 theta(gb,2) = cell_angles(gss).center(gb+1,2)-cell_angles(gss).center(gb,2); %delta column
%                 delta_theta(gb) = atan(theta(gb,1)/theta(gb,2));
%                 delta_phi(gb) = atan(distance(cell_angles(gss).center(gb+1,:),cell_angles(gss).center(gb,:)));
%             elseif cell_angles(gss).cell_slices(1) == cell_angles(gss).cell_slices(1)
%                 theta(gb,:) = [0,0];
%                 delta_theta(gb) = 0;
%                 delta_phi(gb) = 0;
%             else
%                 theta(gb,:) = theta(gb-1,:);
%                 delta_theta(gb) = delta_theta(gb-1);
%                 delta_phi(gb) = delta_phi(gb-1);
%             end
%         end
%         
%         cell_angles(gss).delta_theta = mean(delta_theta);
%         cell_angles(gss).delta_phi = mean(delta_phi);
%     end
%     
%     %%
%     for i = 1:length(cell_angles)-1
%         for j = i+1:length(cell_angles)
%             if abs(cell_angles(i).delta_phi-cell_angles(j).delta_phi) < .1
%                 
%                 first_last = cell_angles(i).center(length(cell_angles(i).cell_slices(2)),:); % Last slice of first cell
%                 last_first = cell_angles(j).center(1,:); % First slice of second cell
%                 slice_diff = cell_angles(j).cell_slices(1)-cell_angles(i).cell_slices(2);
%                 
%                 if slice_diff < 10 && slice_diff > 0
%                     delta_row_i = slice_diff*tan(cell_angles(i).delta_phi)*cos(cell_angles(i).delta_theta);
%                     delta_col_i = slice_diff*tan(cell_angles(i).delta_phi)*sin(cell_angles(i).delta_theta);
%                     delta_row_j = -slice_diff*tan(cell_angles(j).delta_phi)*cos(cell_angles(j).delta_theta);
%                     delta_col_j = -slice_diff*tan(cell_angles(j).delta_phi)*sin(cell_angles(j).delta_theta);
%                     
%                     ix_coord = max([int8(round(first_last(1)+delta_row_i)),1]);
%                     iy_coord = max([int8(round(first_last(2)+delta_col_i)),1]);
%                     jx_coord = max([int8(round(last_first(1)+delta_row_j)),1]);
%                     jy_coord = max([int8(round(last_first(2)+delta_col_j)),1]);
%                     
%                     if total_cells(ix_coord,iy_coord,cell_angles(j).cell_slices(1)) == j && total_cells(jx_coord,jy_coord,cell_angles(i).cell_slices(2)) == i
%                         total_cells(total_cells == j) = i;
%                         delta_row_i2 = tan(cell_angles(i).delta_phi)*cos(cell_angles(i).delta_theta);
%                         delta_col_i2 = tan(cell_angles(i).delta_phi)*sin(cell_angles(i).delta_theta);
%                         im_temp = total_cells(:,:,cell_angles(i).cell_slices(2)) == i;
%                         [xnow,ynow] = ind2sub(size(im_temp),find(im_temp));
%                         
%                         for gn = cell_angles(i).cell_slices(2)+1:cell_angles(j).cell_slices(1)-1
%                             xnow = xnow + delta_row_i2;
%                             xnow(xnow<1) = 1;
%                             ynow = ynow + delta_col_i2;
%                             ynow(ynow<1) = 1;
%                             for ind = 1:length(xnow)
%                                 total_cells(int8(round(xnow(ind))),int8(round(ynow(ind))),gn) = i;
%                             end
%                         end
%                         
%                     end
%                 end
%             end
%         end
%     end
%     
%     
%     
%     for i = 1:length(cell_angles)
%         for g = cell_angles(i).cell_slices(1):cell_angles(i).cell_slices(2)
%             if sum(sum(total_cells(:,:,g)==i)) > 0
%                 last = g;
%             elseif sum(sum(total_cells(:,:,g)==i)) == 0
%                 sli_diff = g - last;
%                 delta_row_i = sli_diff*tan(cell_angles(i).delta_phi)*cos(cell_angles(i).delta_theta);
%                 delta_col_i = sli_diff*tan(cell_angles(i).delta_phi)*sin(cell_angles(i).delta_theta);
%                 im_temp = total_cells(:,:,last) == i;
%                 [xnow,ynow] = ind2sub(size(im_temp),find(im_temp));
%                 xnow = xnow + delta_row_i;
%                 ynow = ynow + delta_col_i;
%                 xnow(xnow<1) = 1;
%                 ynow(ynow<1) = 1;
%                 
%                 for ind = 1:length(xnow)
%                     total_cells(int8(round(xnow(ind))),int8(round(ynow(ind))),g) = i;
%                 end
%             end
%         end
%     end
%     
%     for rk = 1:index
%         slice_sum = 0;
%         cell_slice = [];
%         for g = 1:slices
%             if sum(sum(total_cells(:,:,g)==rk)) > 0
%                 slice_sum = slice_sum + 1;
%                 cell_slice = [cell_slice, g];
%             end
%             
%         end
%         
%         if slice_sum <= 3
%             total_cells(total_cells == rk) = 0;
%         end
%         
%     end
%     
%     
%     indices = unique(total_cells(:));
%     index = 0;
%     
%     for r = 2:length(indices)
%         index = index + 1 ;
%         total_cells(total_cells == indices(r)) = index;
%     end
%     
%     indices = 1:index;
%     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [center, volume] = cellCenter(total_cells,1);
    Single_volume(trial) = length(volume);
    
end