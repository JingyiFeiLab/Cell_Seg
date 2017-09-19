function [slice,stack_o,stack_red,stack_green,stack_blue,stack_back,slices,red_back,green_back,blue_back] = imFormat(filepath,ref_channel,dim,ref_slice,slices2D)

file = dir([filepath '/*.tif']);

red = 1;
green = 2;
blue = 3;
ref = ref_slice;
lp_thresh = .05;


if isempty(findstr(file(1).name,'c1')) && isempty(findstr(file(1).name,'c2')) && isempty(findstr(file(1).name,'c3'))
    
    num = 1:length(file);
    slices = length(num);
    

    if dim == 2
        if ref_slice < 1 + slices2D
            ref = 1+slices2D;
            slice = ref_slice + 1;
        elseif ref_slice > (slices-slices2D)
            ref = slices-slices2D;
            diff_slice = slices - ref_slice;
            slice = ref_slice - 1;
        else
            ref = ref_slice;
            slice = ref_slice - 1;
        end
    end
    I_temp = imread([filepath '/' file(1).name]);
    %I_temp = imread([filepath '.tif']);
    
    channels = size(I_temp,3);

    
    
    if channels == 1
        ref_channel = 1;
        if dim == 3
            stack_o = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_back = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_red = [];
            stack_green = [];
            stack_blue = [];
            red_back = [];
            green_back = [];
            stack_blue = [];
        elseif dim == 2
            stack_o = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_red = [];
            stack_green = [];
            stack_blue = [];
            red_back = [];
            green_back = [];
            stack_blue = [];
        end
    elseif channels == 2
        if dim == 3
            stack_o = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_back = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_red = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_green = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_blue = [];
            red_back = zeros(size(I_temp,1),size(I_temp,2),slices);
            green_back = zeros(size(I_temp,1),size(I_temp,2),slices);
            blue_back = [];
        elseif dim == 2
            stack_o = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_red = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_green = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_blue = [];
            red_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            green_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            blue_back = [];
        end
        
    elseif channels == 3
        if dim == 3
            stack_o = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_back = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_red = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_green = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_blue = zeros(size(I_temp,1),size(I_temp,2),slices);
            red_back = zeros(size(I_temp,1),size(I_temp,2),slices);
            green_back = zeros(size(I_temp,1),size(I_temp,2),slices);
            blue_back = zeros(size(I_temp,1),size(I_temp,2),slices);
        elseif dim == 2
            stack_o = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_red = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_green = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_blue = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            red_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            green_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            blue_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        end
        
    end
    
    if dim == 2
        g = (ref-slices2D):(ref+slices2D);
        slice = find(g==slice);
    elseif dim == 3
        slice = 1:slices;
        g = num;
    end
    
    g_unit = 1;
    for i = g
        I=imread([filepath '/' file(num(i)).name]);
        %I=imread([filepath '.tif']);
        stack_o(:,:,g_unit) = I(:,:,ref_channel);
        if ~isempty(stack_red)
            stack_red(:,:,g_unit) = I(:,:,red);
            red_back(:,:,g_unit) = low_pass(mat2gray(stack_red(:,:,g_unit)),lp_thresh);
        end
        if ~isempty(stack_green)
            stack_green(:,:,g_unit) = I(:,:,green);
            green_back(:,:,g_unit) = low_pass(mat2gray(stack_green(:,:,g_unit)),lp_thresh);
        end
        if ~isempty(stack_blue)
            stack_blue(:,:,g_unit) = I(:,:,blue);
            blue_back(:,:,g_unit) = low_pass(mat2gray(stack_blue(:,:,g_unit)),lp_thresh);
        end
        clear I
        g_unit = g_unit + 1;
    end
    
    pixel_max = max(stack_o(:));
    stack_o = mat2gray(stack_o, [0,double(pixel_max)]);

else
    
    
    I_temp = imread([filepath '/' file(1).name]);
    number = length(file);
    slices = number/3;
    num1 = 1:3:number;
    num2 = 2:3:number;
    num3 = 3:3:number;
    
    if ref_channel == 1
        sharks = num1;
    elseif ref_channel == 2
        sharks = num2;
    elseif ref_channel == 3
        sharks = num3;
    end
    
    if dim == 2
        if ref_slice < 1 + slices2D
            ref = 1+slices2D;
            slice = ref_slice + 1;
        elseif ref_slice > (slices-slices2D)
            ref = slices-slices2D;
            diff_slice = slices - ref_slice;
            slice = ref_slice - 1;
        else
            ref = ref_slice;
            slice = ref_slice - 1;
        end
        
    
        yup = (ref-slices2D):(ref+slices2D);
        stack_o = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        stack_red = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        stack_green = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        stack_blue = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        stack_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        red_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        green_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        blue_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        slice = find(yup==slice);
    elseif dim == 3
        slice = 1:slices;
        yup = slice;
        stack_o = zeros(size(I_temp,1),size(I_temp,2),slices);
        stack_red = zeros(size(I_temp,1),size(I_temp,2),slices);
        stack_green = zeros(size(I_temp,1),size(I_temp,2),slices);
        stack_blue = zeros(size(I_temp,1),size(I_temp,2),slices);
        stack_back = zeros(size(I_temp,1),size(I_temp,2),slices);
        red_back = zeros(size(I_temp,1),size(I_temp,2),slices);
        green_back = zeros(size(I_temp,1),size(I_temp,2),slices);
        blue_back = zeros(size(I_temp,1),size(I_temp,2),slices);
    end
    
    g_unit = 1;
    for g = yup
        stack_o(:,:,g_unit)=imread([filepath '/' file(sharks(g)).name]);
        stack_red(:,:,g_unit) = imread([filepath '/' file(num1(g)).name]);
        stack_green(:,:,g_unit) = imread([filepath '/' file(num2(g)).name]);
        stack_blue(:,:,g_unit) = imread([filepath '/' file(num3(g)).name]);
        red_back(:,:,g_unit) = low_pass(mat2gray(stack_red(:,:,g_unit)),lp_thresh);
        green_back(:,:,g_unit) = low_pass(mat2gray(stack_green(:,:,g_unit)),lp_thresh);
        blue_back(:,:,g_unit) = low_pass(mat2gray(stack_blue(:,:,g_unit)),lp_thresh);
        
        g_unit = g_unit + 1;
    end
    
    pixel_max = max(stack_o(:));
    stack_o = mat2gray(stack_o, [0, double(pixel_max)]);
    
end

red_back = red_back * double(max(stack_red(:)));
green_back = green_back * double(max(stack_green(:)));
blue_back = blue_back * double(max(stack_blue(:)));


end