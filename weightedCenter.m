function w_center = weightedCenter(I,objects,id)

[row,col] = size(objects);
object_i = objects == id;
pixels = sub2ind(size(object_i),find(object_i));
w_center = [0,0];
mass = sum(sum(object_i.*I));

for i = 1:length(pixels)
    [x,y] = ind2sub(size(object_i),pixels(i));
    w_center(1,1) = w_center(1,1) + x*I(x,y);
    w_center(1,2) = w_center(1,2) + y*I(x,y);
end

w_center = w_center/mass;



end