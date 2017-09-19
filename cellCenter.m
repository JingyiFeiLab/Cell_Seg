function [center,pixels] = cellCenter(objects, id)

if ndims(objects) == 2
    [row,col] = size(objects);
    object_i = objects == id;
    pixels = sub2ind(size(object_i),find(object_i));
    center = [0,0];
    
    for i = 1:length(pixels)
        [x,y] = ind2sub(size(object_i),pixels(i));
        center(1,1) = center(1,1) + x;
        center(1,2) = center(1,2) + y;
    end
    
    center = floor(center)/length(pixels);

elseif ndims(objects) == 3;
    [row,col, v] = size(objects);
    object_i = objects == id;
    pixels = sub2ind(size(object_i),find(object_i));
    center = [0,0,0];
    
    for i = 1:length(pixels)
        [x,y,z] = ind2sub(size(object_i),pixels(i));
        center(1,1) = center(1,1) + x;
        center(1,2) = center(1,2) + y;
        center(1,3) = center(1,3) + z;
    end
    
    center = floor(center)/length(pixels);

end