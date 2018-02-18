function [slice,stack_o,stack_one,stack_two,stack_three,stack_back,slices] = imFormat(filepath,ref_channel,dim,ref_slice,slices2D,channels)

file = dir([filepath '/*.tif']);
%file = strcat([filepath '.tif']);

one = 1;
two = 2;
three = 3;
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
            stack_o = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_one = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_two = [];
            stack_three = [];
            
        elseif dim == 2
            stack_o = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_one = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));;
            stack_two = [];
            stack_three = [];
            
        end
    elseif channels == 2
        if dim == 3
            stack_o = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_back = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_one = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_two = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_three = [];
            
        elseif dim == 2
            stack_o = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_one = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_two = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_three = [];
            
        end
        
    elseif channels == 3
        if dim == 3
            stack_o = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_back = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_one = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_two = zeros(size(I_temp,1),size(I_temp,2),slices);
            stack_three = zeros(size(I_temp,1),size(I_temp,2),slices);
            
        elseif dim == 2
            stack_o = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_one = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_two = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            stack_three = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
            
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
        if ~isempty(stack_one)
            stack_one(:,:,g_unit) = I(:,:,one);
            one_back(:,:,g_unit) = low_pass(mat2gray(stack_one(:,:,g_unit)),lp_thresh);
        end
        if ~isempty(stack_two)
            stack_two(:,:,g_unit) = I(:,:,two);
            two_back(:,:,g_unit) = low_pass(mat2gray(stack_two(:,:,g_unit)),lp_thresh);
        end
        if ~isempty(stack_three)
            stack_three(:,:,g_unit) = I(:,:,three);
            three_back(:,:,g_unit) = low_pass(mat2gray(stack_three(:,:,g_unit)),lp_thresh);
        end
        clear I
        g_unit = g_unit + 1;
    end
    
    pixel_max = max(stack_o(:));
    stack_o = mat2gray(stack_o, [0,double(pixel_max)]);

else
    
    
    I_temp = imread([filepath '/' file(1).name]);
    number = length(file);
    slices = number/channels;
    num1 = 1:channels:number;
    num2 = 2:channels:number;
    num3 = 3:channels:number;
    
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
        stack_one = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        stack_two = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        stack_three = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        stack_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        one_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        two_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        three_back = zeros(size(I_temp,1),size(I_temp,2),(2*slices2D+1));
        slice = find(yup==slice);
    elseif dim == 3
        slice = 1:slices;
        yup = slice;
        stack_o = zeros(size(I_temp,1),size(I_temp,2),slices);
        stack_one = zeros(size(I_temp,1),size(I_temp,2),slices);
        stack_two = zeros(size(I_temp,1),size(I_temp,2),slices);
        stack_three = zeros(size(I_temp,1),size(I_temp,2),slices);
        stack_back = zeros(size(I_temp,1),size(I_temp,2),slices);
        one_back = zeros(size(I_temp,1),size(I_temp,2),slices);
        two_back = zeros(size(I_temp,1),size(I_temp,2),slices);
        three_back = zeros(size(I_temp,1),size(I_temp,2),slices);
    end
    
    g_unit = 1;
    for g = yup
        stack_o(:,:,g_unit)=imread([filepath '/' file(sharks(g)).name]);
        stack_one(:,:,g_unit) = imread([filepath '/' file(num1(g)).name]);
        stack_two(:,:,g_unit) = imread([filepath '/' file(num2(g)).name]);
        stack_three(:,:,g_unit) = imread([filepath '/' file(num3(g)).name]);
        one_back(:,:,g_unit) = low_pass(mat2gray(stack_one(:,:,g_unit)),lp_thresh);
        two_back(:,:,g_unit) = low_pass(mat2gray(stack_two(:,:,g_unit)),lp_thresh);
        three_back(:,:,g_unit) = low_pass(mat2gray(stack_three(:,:,g_unit)),lp_thresh);
        
        g_unit = g_unit + 1;
    end
    
    pixel_max = max(stack_o(:));
    stack_o = mat2gray(stack_o, [0, double(pixel_max)]);
    
end



end