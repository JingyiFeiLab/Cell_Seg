%Low Pass filter for cell image stacks
%Matt Reyer 5/23/16

function im_lp = low_pass(image, fc) % image = single layer original image , fc = cutoff frequency [0,1] 
    [ir,ic] = size(image);
    hr = (ir-1)/2;
    hc = (ic-1)/2;
    [x,y] = meshgrid(-hc:hc, -hr:hr);
    
    mg = sqrt((x/hc).^2 + (y/hr).^2);
    Ip = double(mg <= fc);
    
    IM = fftshift(fft2(double(image)));
    IP = zeros(size(IM));
    IP(:,:) = IM(:,:).*Ip;
    
    im_lp = abs(ifft2(ifftshift(IP)));
end

    
 
    
    



































