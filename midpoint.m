function midpoint = midpoint(P1,P2)

% Input : 2 Points (x,y,z)
% Output : Euclidean Distance

if length(P1) == 2 % 2D
    x1 = P1(2);
    x2 = P2(2);
    y1 = P1(1);
    y2 = P2(1);
    
    midpoint = [(x1+x2)/2,(y1+y2)/2];
    
elseif length(P1) == 3 % 3D
    
    x1 = P1(2);
    x2 = P2(2);
    
    y1 = P1(1);
    y2 = P2(1);
    
    z1 = P1(3);
    z2 = P2(3);
    
    midpoint = [(x1+x2)/2,(y1+y2)/2,(z1+z2)/2];
    
    
end

end