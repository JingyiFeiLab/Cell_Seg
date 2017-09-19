function im = echbreak(bwim);

for i = 2:length(bwim)-1
        for j = 2:length(bwim)-1
            x_neighbors = bwim(i-1,j) + bwim(i+1,j);
            y_neighbors = bwim(i, j-1) + bwim(i, j+1);
            if x_neighbors == 0 || y_neighbors == 0
                bwim(i,j) = 0;
            end
        end
end
im = bwim;
end
