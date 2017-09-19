function ellipticity = cellEllipse(objects,id, varargin)

ellipticity = zeros(1,4);

ellipse = objects == id;
[X Y] = ind2sub(size(ellipse),find(ellipse));
E = fit_ellipse(X,Y);
[ellipticity(1,1), ellipticity(1,2), ellipticity(1,3), ellipticity(1,4)] = find_ellipse(E);


if ellipticity(1,1) == 0
    return
end


ellipticity;


end