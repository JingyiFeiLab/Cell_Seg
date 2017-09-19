function error = ellipseTest(edge1,edge2,area,pix_size)

% INPUT : edge1 = edge pixels of interest (in cell array)
%         edge2 = edge pixels of test ellipse (also in cell array)
%
% OUTPUT : error = sum of squared distances between edge1 pixels and
%                  nearest edge2 pixel

if isempty(edge1) || isempty(edge2)
    error = 100;
    return
end

pix = length(edge1{1,1});
pix2 = length(edge2{1,1});

area = area/(pix_size^2);

differences1 = zeros(pix,1);
differences2 = zeros(pix2,1);

for i = 1:pix
    distances = zeros(pix2,1);
    for j = 1:pix2
        distances(j) = distance(edge1{1,1}(i,:),edge2{1,1}(j,:));
    end
    differences1(i) = min(distances)^2;
end

for k = 1:pix2
    distances2 = zeros(pix,1);
    for l = 1:pix
        distances2(l) = distance(edge1{1,1}(l,:),edge2{1,1}(k,:));
    end
    differences2(k) = min(distances2)^2;
end


error1 = sum(differences1)/pix2;
error2 = sum(differences2)/pix;

error = max(error1,error2);

end