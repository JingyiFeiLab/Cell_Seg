function [ellipse1,test] = ellipseError(objects, id)
% objects = Labeled binary image
% id = cell/non-cell in question

ellipse = objects == id;
ellipse1 = bwboundaries(ellipse);
%[X1,Y1] = ind2sub(size(ellipse1),find(ellipse1));
[X,Y] = ind2sub(size(ellipse),find(ellipse));
E = fit_ellipse(X,Y);

if isfield(E,'angleToX') == 0
    ellipse1 = {};
    test = {};
    
elseif E.long_axis > 100
    ellipse1 = {};
    test = {};
else
    if E.angleToX == E.angleFromX
        [x_test,y_test] = calcEllipse(E.X0_in,E.Y0_in,E.a,E.b,E.angleToX,360);
    elseif E.angleToX ~= E.angleFromX
        [x_test,y_test] = calcEllipse(E.X0_in,E.Y0_in,E.a,E.b,180-E.angleToX,360);   
    end

    test = zeros(size(objects));
    
    if sum(x_test==0)>0 || sum(y_test==0) > 0 
        ellipse1 = {};
        test = {};
        return
    end
    
    for i = 1:length(x_test)
        test(x_test(i),y_test(i)) = 1;
    end
    
    test = imfill(test,'holes');
    test = bwboundaries(test);
end


end