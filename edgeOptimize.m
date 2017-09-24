function [edges,new_bounds,con_peaks] = edgeOptimize(BW,id)

%BW_label = bwlabel(BW,4);
%objects = max(BW_label(:));


new_bounds = cell(1,1);
con_peaks = zeros(1,1);
edges = zeros(size(BW));

for i = id
%for i = 3
    %i
    w = 13;
    peaks = 0;
    round = 0;
    
    while peaks < 2 && round < 2
        round = round + 1;
        %BWi = BW_label == i;
        BWi = BW == i;
        BWi_bounds = bwperim(BWi);
        %BWi_bounds = edgeBreak(BWi_bounds);
        %bounds = bwboundaries(BW_label==i);
        bounds = bwboundaries(BW==i);
        if length(bounds{1,1}) < w-4
            continue
        end
        w = w - 4; %Savitzky-Golay Filter Frame Width
        p = 2;  %Savitzky-Golay Filter Polynomial Order
        %     new_bounds{i,1} = zeros(length(bounds{1,1}),7);
        %     new_bounds{i,1}(:,1) = bounds{1,1}(:,1);
        %     new_bounds{i,1}(:,2) = bounds{1,1}(:,2);
        new_bounds{1,1}(:,1) = sgolayfilt(bounds{1,1}(:,1),p,w);
        new_bounds{1,1}(:,2) = sgolayfilt(bounds{1,1}(:,2),p,w);
        
        
        l = length(new_bounds{1,1}(:,1));
        for j = 1:l
            p1 = new_bounds{1,1}(j,:);
            if j == 1 || j == 2
                p_minus1 = new_bounds{1,1}(l-2,:);
            else
                p_minus1 = new_bounds{1,1}(mod(j-1,l),:);
            end
            
            if j == l-1 || j == l-2
                p_plus1 = new_bounds{1,1}(1,:);
            else
                p_plus1 = new_bounds{1,1}(mod(j+1,l),:);
            end
            
            if j == 1 || j == 2
                p_minus2 = new_bounds{1,1}(l-2,:);
            else
                p_minus2 = new_bounds{1,1}(mod(j-2,l),:);
            end
            
            if j == l-2
                p_plus2 = new_bounds{1,1}(1,:);
            else
                p_plus2 = new_bounds{1,1}(mod(j+2,l),:);
            end
            
            new_bounds{1,1}(j,3) = calcNormal(p_minus1([1,2]),p_plus1([1,2]));
        end
        
        for j = 1:l
        %for j = 36
            p1 = new_bounds{1,1}(j,:);
            if j == 1 || j == 2
                p_minus1 = new_bounds{1,1}(l-2,:);
            else
                p_minus1 = new_bounds{1,1}(j-1,:);
            end
            
            
            if p1(3) == 0
                [~,I] = min([abs(2*pi-p_minus1(3)),abs(0-p_minus1(3))]);
                if I == 1
                    p1(3) = 2*pi;
                else
                    p1(3) = 0;
                end
            end
            
            if p_minus1(3) == 0
                [~,I] = min([abs(p1(3) -2*pi),abs(p1(3)-0)]);
                if I == 1
                    p_minus1(3) = 2*pi;
                else
                    p_minus1(3) = 0;
                end
            end
            
            if (p1(3) < 1 && p1(3) > 0 && p_minus1(3) < 2*pi && p_minus1(3) > 2*pi - 1) 
                new_bounds{1,1}(j,4) = mod(p1(3)-p_minus1(3),2*pi);
            elseif (p_minus1(3) < 1 && p_minus1(3) > 0 && p1(3) < 2*pi && p1(3) > 2*pi - 1)
                new_bounds{1,1}(j,4) = -1*mod(p_minus1(3)-p1(3),2*pi);
            else
                new_bounds{1,1}(j,4) = mod((p1(3) - p_minus1(3))/sign(p1(3)-p_minus1(3)),2*pi)/sign(p1(3)-p_minus1(3));
            end
            new_bounds{1,1}(isnan(new_bounds{1,1}(:,4)),4) = 0;
        end
        
        new_bounds{1,1}(:,4) = smooth(new_bounds{1,1}(:,4),3);
        cave_points = new_bounds{1,1}(:,4);
        cave_points(cave_points<=0.1) = 0;
        [peaks,~] = findpeaks(cave_points);
        peaks = length(peaks);
        
    end
    if exist('cave_points') == 1
        peaks1 = findpeaks(cave_points);
        con_peaks(1,1) = sum(peaks1>=.3);
    else
        con_peaks = 0;
    end
end



end