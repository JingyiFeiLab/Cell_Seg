function [split_im] = concaveSplit(BW, id, pix_size,ee_thresh);

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
    cave_points(cave_points<=0.1) = 0;
    [peaks,locs] = findpeaks(cave_points);
    max_cave = [peaks, locs];
    

    
    if length(max_cave(:,1)) < 1 | size(max_cave) == [0,0] | size(max_cave) == [1,2];
        split_im = object;
        break
    end
    [~,order] = sort(max_cave(:,1),'descend');
    max_cave = max_cave(order,:);
    num_cave = length(max_cave);
    split_lines = {};
    line_index = 1;
    
    for j = 1:num_cave-1
        for k = j+1:num_cave
            
            if num_cave == 1
                break
            end
            
            x_concave = [reg_bounds{1,1}(max_cave(j,2),2),reg_bounds{1,1}(max_cave(k,2),2)];
            y_concave = [reg_bounds{1,1}(max_cave(j,2),1),reg_bounds{1,1}(max_cave(k,2),1)];
            
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
            split_lines{line_index,1} = line;
            line_index = line_index + 1;
            
        end
        
        clear x y line x_concave y_concave slope
        
    end
    
    split_lines = split_lines(~cellfun('isempty',split_lines));
    split_num = length(split_lines) + 1;
    
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
        split_lines{line_index,1} = line;
    end
    
    split_lines = split_lines(~cellfun('isempty',split_lines));
    object_composite = zeros(size(object));
    
    while 0 < 1
        error_array = [0,0,100];
        final_error_array = [0,100];
        [ellipse_old,test_old] = ellipseError(object,1);
        if isempty(ellipse_old) == 1 || isempty(test_old) == 1
            error_old = ee_thresh+1;
        else
            error_old = ellipseTest(ellipse_old,test_old,cellArea(object,1),pix_size);
        end
        
        if error_old < ee_thresh
            object_composite = object_composite + object;
            break
        end
        
        for m = 1:length(split_lines)
            
            object1 = object;
            object1(sub2ind(size(object),split_lines{m}(:,1),split_lines{m}(:,2))) = 0;
            object1 = bwareaopen(object1,20,4);
            new_objects = bwlabel(object1,4);
            num_new = max(new_objects(:));
            e_error = zeros(1,num_new);
            
            
            for n = 1:num_new
                [ellipse1, test] = ellipseError(new_objects,n);
                if isempty(ellipse1) == 1 || isempty(test) == 1
                    e_error(n) = ee_thresh+1;
                else
                    e_error(n) = ellipseTest(ellipse1,test,cellArea(new_objects,n),pix_size);
                end
            end
            
            [min_error,mindex] = min(e_error);
            
            if length(e_error) == 2 && e_error(1) < ee_thresh && e_error(2)< ee_thresh && sum(e_error) < final_error_array(2)
                final_error_array(1) = m;
                final_error_array(2) = sum(e_error);
            end
            
            if min_error < error_array(3)
                error_array(1) = m;
                error_array(2) = mindex;
                error_array(3) = min_error;
            end
            
        end
        
        if length(e_error) == 2 && e_error(1) < ee_thresh && e_error(2)< ee_thresh && sum(e_error) < final_error_array(2)
            object_composite = object_composite_old;
            if isempty(split_lines)
                break
            end
            object_new = object_new_old;
            object_new(sub2ind(size(object),split_lines{error_array(1)}(:,1),split_lines{error_array(1)}(:,2))) = 0;
            object_new2 = bwlabel(object_new,4);
            object_new_new = (object_new2 == error_array(2));
            object_old = (object_new ~= object_new_new);
            object = object_old;
            object_composite_old = object_composite;
            object_composite = object_composite + object_new_new;
            split_lines{error_array(1)} = [];
            split_lines = split_lines(~cellfun('isempty',split_lines));
            break
        elseif error_array(3) < ee_thresh && error_array(3) < error_old
            
            object_new = object;
            object_new_old = object_new;
            object_new(sub2ind(size(object),split_lines{error_array(1)}(:,1),split_lines{error_array(1)}(:,2))) = 0;
            object_new2 = bwlabel(object_new,4);
            object_new_new = (object_new2 == error_array(2));
            object_old = (object_new ~= object_new_new);
            object = object_old;
            object_composite_old = object_composite;
            object_composite = object_composite + object_new_new;
            split_lines{error_array(1)} = [];
            split_lines = split_lines(~cellfun('isempty',split_lines));
        else
            
            break
        end
        
        
    end
    
    split_im = object_composite;

end