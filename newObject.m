function [new_objects,new_mask] = newObject(object,mask,new_borders);

max_id = max(max(object));
new_object = zeros(size(object,1),size(object,2));

for i = 1:length(new_borders)
    for j = 1:length(new_borders(i).New_Borders)
        new_object(new_borders(i).New_Borders(j,2),new_borders(i).New_Borders(j,1)) = max_id+1;
    end
    new_object = imfill(new_object,'holes');
    max_id = max_id + 1;
end

object = object + new_object;
new_mask = mask + im2bw(new_object,.99);

new_objects = object;

end