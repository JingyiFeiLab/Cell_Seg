function [split_im, split_lines] = concaveSplit(BW, id, pix_size,ee_thresh);

% INPUT : BW = Labeled BW image
% OUTPUT : split_im = Image with Non-Single Cells split at concave points



[edge,bounds,~] = edgeOptimize(BW,id);
xdim = size(BW,1);
ydim = size(BW,2);

%objects = max(BW(:));


%for i = 1:objects
for i = id     
    
    object = BW == i;
    reg_bounds = bwboundaries(object);
    cave_points = bounds{1,1}(:,4);
    cave_points(cave_points<=0) = 0;
    [peaks,locs] = findpeaks(cave_points);
    max_cave = [peaks, locs];
%     for i = 1:length(max_cave)
%         max_cave(i,2) = max_cave(i,2)-1;
%         if max_cave(i,2) == 0
%             max_cave(i,2) = length(bounds{i,1});
%         end
%     end
    
    if length(max_cave(:,1)) < 2 | size(max_cave) == [0,0];
        break
    end
    [~,order] = sort(max_cave(:,1),'descend');
    max_cave = max_cave(order,:);
    num_cave = length(max_cave);
    split_lines = cell(num_cave,num_cave);
    
    
    
    if max(max_cave(:,1)) > .35 && sum(max_cave(:,1)>.15) < 2 || max(max_cave(:,1)) > .5
        
        x_concave = reg_bounds{1,1}(max_cave(1,2),2);
        y_concave = reg_bounds{1,1}(max_cave(1,2),1);
        
        pointx = x_concave;
        pointy = y_concave;
        
        slope = -1*(1/tan(bounds{1,1}(max_cave(1,2),3)));
        
        if (abs(slope) <= .1) && (bounds{1,1}(max_cave(1,2),3) >= 1.5 * (pi-.2) || bounds{1,1}(max_cave(1,2),3) <= 1.5 * (pi+.2))
            while object(pointy,pointx) == 1
                pointy = pointy - 1;
            end
        elseif (abs(slope) <= .1) && (bounds{1,1}(max_cave(1,2),3) >= .5 * (pi-.2) || bounds{1,1}(max_cave(1,2),3) <= .5 * (pi+.2))
            while object(pointy,pointx) == 1
                pointy = pointy + 1;
            end 
        elseif (abs(slope) >= 30) && (bounds{1,1}(max_cave(1,2),3) >= (pi-.2) || bounds{1,1}(max_cave(1,2),3) <= (pi+.2))
               while object(pointy,pointx) == 1
                pointx = pointx + 1;
               end 
            
        elseif (abs(slope) >= 30) && (bounds{1,1}(max_cave(1,2),3) >= 2*(pi-.2) || bounds{1,1}(max_cave(1,2),3) <= 2)  
            while object(pointy,pointx) == 1
                pointx = pointx - 1;
            end
             
        elseif bounds{1,1}(max_cave(1,2),3) > pi && bounds{1,1}(max_cave(1,2),3) < 1.5*pi
            while object(pointy,pointx) == 1
                
                pointx = int32(round(pointx + (1/abs(slope))));
                pointy = pointy - 1;
                
                if pointx > xdim
                    pointx = xdim;
                end
                
            end
            
        elseif bounds{1,1}(max_cave(1,2),3) > 0 && bounds{1,1}(max_cave(1,2),3) < .5*pi
            while object(pointy,pointx) == 1
                
                pointx = int32(round(pointx - (1/abs(slope))));
                pointy = pointy + 1;
                
                if pointx < 1
                    pointx = 1;
                end
                
            end
        
        elseif bounds{1,1}(max_cave(1,2),3) > .5*pi && bounds{1,1}(max_cave(1,2),3) < pi
            while object(pointy,pointx) == 1
                
                pointx = int32(round(pointx + (1/abs(slope))));
                pointy = pointy + 1;
                
                if pointx > xdim
                    pointx = xdim;
                end
                
            end    
         
        elseif bounds{1,1}(max_cave(1,2),3) > 1.5*pi && bounds{1,1}(max_cave(1,2),3) < 2*pi
            while object(pointy,pointx) == 1
                
                pointx = int32(round(pointx - (1/abs(slope))));
                pointy = pointy - 1;
                
                if pointx < 1
                    pointx = 1;
                end
                
            end   
               
        end
        
        slope = (pointy-y_concave)/(pointx-x_concave);
        
        if abs(slope) >= 30
            y = min(y_concave,pointy):max(y_concave,pointy);
            x = repmat(x_concave,1,length(y));
            y = y';
            x = x';
            
            line = [y,x];
            
            
        elseif abs(slope) <= .1 
            x = min(x_concave,pointx):max(x_concave,pointx);
            y = repmat(y_concave,1,length(x));
            x = x';
            y = y';
            line = [y,x];
            
        else
            x_dist = abs(x_concave-pointx);
            y_dist = abs(y_concave-pointy);
            if x_dist >= y_dist
                x = linspace(min(x_concave,double(pointx)),max(x_concave,double(pointx)),distance([x_concave y_concave], [double(pointx),double(pointy)])+3) ;
                x = x';
                for p = 1:length(x(:,1))
                    y(p,1) = slope*(x(p,1)-x_concave) + y_concave ;
                end
            elseif y_dist > x_dist
                
                y = linspace(min(y_concave,double(pointy)),max(y_concave,double(pointy)),distance([x_concave y_concave], [double(pointx),double(pointy)])+3) ;
                y = y';
                for p = 1:length(y(:,1))
                    x(p,1) = (y(p,1) - y_concave)/slope + x_concave;
                end
            end
            
            x = floor(x);
            y = floor(y);
            
            line = [y,x];
            
        end
        
        line(line<1) = 1;
        line(line>xdim) = xdim;
        line(line>ydim) = ydim;
        split_lines{1,1} = line;
    
    else
        
        for j = 1:num_cave-1
            for k = j+1:num_cave
                
                edge_pic1 = zeros(length(reg_bounds{1,1}(:,1)),1);
                edge_pic2 = zeros(length(reg_bounds{1,1}(:,1)),1);
                
                for jd = 1:length(reg_bounds{1,1}(:,1))
                    edge_pic1(jd) = distance(reg_bounds{1,1}(jd,:),bounds{1,1}(max_cave(j,2),:));
                end
                
                [~,e1] = min(edge_pic1);
                
                for kd = 1:length(reg_bounds{1,1}(:,1))
                    edge_pic2(kd) = distance(reg_bounds{1,1}(kd,:),bounds{1,1}(max_cave(k,2),:));
                end
                
                [~,e2] = min(edge_pic2);
                
                x_concave = [reg_bounds{1,1}(e1,2),reg_bounds{1,1}(e2,2)];
                y_concave = [reg_bounds{1,1}(e1,1),reg_bounds{1,1}(e2,1)];
                
                if x_concave(1) == x_concave(2) && y_concave(1) == y_concave(2)
                    continue
                end
                
                slope = zeros(2,2);
                for e =1:length(x_concave)
                    for f= 1:length(y_concave)
                        slope(e,f) = (y_concave(e)-y_concave(f))/(x_concave(e)-x_concave(f));
                    end
                end
                
                %Connect Concave Points
                for h = 1
                    for o= 2
                        if slope(h,o) == -Inf || slope(h,o) == Inf
                            y = min(y_concave(h),y_concave(o)):max(y_concave(h),y_concave(o));
                            x = repmat(x_concave(h),1,length(y));
                            y = y';
                            x = x';
                            
                            line = [y,x];
                            
                            
                        elseif slope(h,o) == 0
                            x = min(x_concave(h),x_concave(o)):max(x_concave(h),x_concave(o));
                            y = repmat(y_concave(h),1,length(x));
                            x = x';
                            y = y';
                            line = [y,x];
                            
                        else
                            x_dist = abs(x_concave(h)-x_concave(o));
                            y_dist = abs(y_concave(h)-y_concave(o));
                            if x_dist >= y_dist
                                x = linspace(min(x_concave(h),x_concave(o)),max(x_concave(h),x_concave(o)),distance([x_concave(h) y_concave(h)], [x_concave(o),y_concave(o)])+1) ;
                                x = x';
                                for p = 1:length(x(:,1))
                                    y(p,1) = slope(h,o)*(x(p,1)-x_concave(h)) + y_concave(h) ;
                                end
                            elseif y_dist > x_dist
                                y = linspace(min(y_concave(h),y_concave(o)),max(y_concave(h),y_concave(o)),distance([x_concave(h) y_concave(h)], [x_concave(o),y_concave(o)])+1) ;
                                y = y';
                                for p = 1:length(y(:,1))
                                    x(p,1) = (y(p,1) - y_concave(h))/slope(h,o) + x_concave(h);
                                end
                            end
                            
                            x = floor(x);
                            y = floor(y);
                            
                            line = [y,x];
                            
                            
                        end
                        
                        y = [];
                        x= [];
                        
                    end
                    
                end
                
                line(line<1) = 1;
                line(line>xdim) = xdim;
                line(line>ydim) = ydim;
                split_lines{j,k} = line;
                
            end
            
            clear x y line x_concave y_concave slope
            
        end
    end
    
    split_lines = split_lines(~cellfun('isempty',split_lines));
    
    for m = 1:length(split_lines)
        
        object1 = object;
        object1(sub2ind(size(object),split_lines{m}(:,1),split_lines{m}(:,2))) = 0;
        object1 = bwareaopen(object1,20,4);
        new_objects = bwlabel(object1,4);
        num_new = max(new_objects(:));
        prob = zeros(num_new,1);
        e_error = zeros(num_new,1);
        single_cells = [];
        
        
        for n = 1:num_new
            
            ellipticity = cellEllipse(new_objects,n);
            area = cellArea(new_objects,n);
            [ellipse1, test] = ellipseError(new_objects,n);
            if isempty(ellipse1) == 1 || isempty(test) == 1
                e_error(n) = ee_thresh+1;
            else
                e_error(n) = ellipseTest(ellipse1,test,area,pix_size);
            end
            %prob(n) = cell_prob(ellipticity(1,1),area,ellip_compare,area_compare,density);
            
            if e_error(n) < ee_thresh %&& prob(n) > thresh
                
                single_cells = [single_cells 1];
                break
            
            end
            
            if sum(single_cells) >= 1
                break
            end
            
        end
        
        if sum(single_cells) >= 1
            break 
        else 
            if m > 1
                for md = 1:m-1
                    object1 = object;
                    
                    object1(sub2ind(size(object),split_lines{m}(:,1),split_lines{m}(:,2))) = 0;
                    object1(sub2ind(size(object),split_lines{md}(:,1),split_lines{md}(:,2))) = 0;
                    object1 = bwareaopen(object1,20,4);
                    new_objects = bwlabel(object1,4);
                    num_new = max(new_objects(:));
                    %prob = zeros(num_new,1);
                    e_error = zeros(num_new,1);
                    single_cells = [];
                    
                    
                    for nd = 1:num_new
                        ellipticity = cellEllipse(new_objects,nd);
                        area = cellArea(new_objects,nd,pix_size);
                        [ellipse1, test] = ellipseError(new_objects,nd);
                        if isempty(ellipse1) == 1 || isempty(test) == 1
                            e_error(nd) = 2;
                        else
                            e_error(nd) = ellipseTest(ellipse1,test,area,pix_size);
                        end
                        
                        %prob(nd) = cell_prob(ellipticity(1,1),area,ellip_compare,area_compare,density);
                        
                        if e_error(nd) < ee_thresh %&& prob(nd) > thresh
                            
                            single_cells = [single_cells 1];
                            break
                            
                        end
                        
                        if sum(single_cells) >= 1
                            break
                        end
                        
                    end
                    
                    if sum(single_cells) >=1
                        break 
                    end
                    
                end
                
                if sum(single_cells) >= 1
                    break
                end
            end
        end
        
        if sum(single_cells) >= 1
            break 
        end
    
        end
 
end

if (exist('single_cells') == 1 && sum(single_cells) >=1) 
    split_im = object1;
else 
    split_im = object;
    split_lines = {};
end
    

end