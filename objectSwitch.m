function new_mask = objectSwitch(objects11,mask01,objects01,mask12,objects12)

% Objects11 = post-manualCorrect Objects (i.e. newfig , part3)
% mask01 = pre-manualCorrect mask (i.e. part2(g).Mask)
% Objects01 = pre-manualCorrect (i.e. part2(g).Objects)
% mask12 = part2(g+1).Mask
% Objects12 = part2(g+1).Objects

%   MANUAL CORRECTION
% If cell x, frame g switched from negative to positive --> set cell x,
% frame g + 1 as positive
%
% If cell x, frame g switched from positive to negative --> set cell x,
% frame g + 1 as negative

max_id = max(max(objects11));
max_id2 = max(max(objects01));

for i = 1:max_id;
    
    % Negative to Positive
    if sum(sum((objects11 == i).*mask01)) == 0 && sum(sum((objects11 == i).*objects01)) ~= 0
        id = max(max(im2bw(objects11 == i,.99).*objects12));
        if sum(sum((objects12 == id).*mask12)) == 0
            mask12 = mask12 + (objects12 == id);
        end
    end
    
end

for j = 1:max_id2
    % Positive to Negative
    if sum(sum((objects01 == j).*objects11)) == 0
        id = max(max(im2bw(objects01 == j).*objects12));
        mask12 = mask12 - (objects12 == id);
    end


new_mask = mask12;

end


