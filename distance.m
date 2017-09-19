function distance = distance(P1, P2)
% Input : 2 Points (x,y)
% Output : Euclidean Distance


if isa(P1,'double') == 0
    P1 = double(P1);
end

if isa(P2,'double') == 0
    P2 = double(P2);
end

if length(P1) == 2 % 2D
    x1 = P1(1);
    x2 = P2(1);
    y1 = P1(2);
    y2 = P2(2);
    
    distance = sqrt((x1-x2)^2 + (y1-y2)^2);

elseif length(P1) == 3 % 3D
    
    x1 = P1(1);
    x2 = P2(1);
    
    y1 = P1(2);
    y2 = P2(2);
    
    z1 = P1(3);
    z2 = P2(3);
    
    distance = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2);
    
    
end

end

