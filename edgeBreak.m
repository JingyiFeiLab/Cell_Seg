function [edges,BW] = edgeBreak(BW)

[row, col] = size(BW);

edges = bwperim(BW);

for i = 2:row-1
    for j = 2:col-1
        if edges(i-1,j)+edges(i+1,j)+edges(i,j-1)+edges(i,j+1) >= 3
            edges(i,j) = 0;
            BW(i,j) = 0;
        end
    end
end

edges;
BW = bwareaopen(BW,20);

end