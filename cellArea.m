function area = cellArea(objects, id, varargin)

% varargin = Pixel Area =  length/pixel

area = sum(objects(:) == id)/1.17; %Constant to account for edge pixels

if nargin == 3
    pix_area = varargin{1};
    pix_area = pix_area^2;
    area = pix_area * area;
end

area;



end