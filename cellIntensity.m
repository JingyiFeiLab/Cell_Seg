function intensity = cellIntensity(objects,I,id)

%[r,c,v] = size(objects);
se = zeros(3,3,3);
se(:,:,2) = [1,1,1;1,1,1;1,1,1];
object_i = objects == id;
object_i = imdilate(object_i,se);
intensity = sum(sum(sum(object_i.*I)));

end





