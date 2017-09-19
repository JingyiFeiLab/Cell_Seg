function line = makeLine(xv,yv)

x_concave = double(xv);
y_concave = double(yv);



slope = zeros(2,2);
for e =1:length(x_concave)
    for f= 1:length(y_concave)
        slope(e,f) = (y_concave(e)-y_concave(f))/(x_concave(e)-x_concave(f));
    end
end

%Connect Concave Points
for h = 1
    for o= 2
        if slope(h,o) == -Inf || slope(h,o) == Inf
            y = min(y_concave(h),y_concave(o)):max(y_concave(h),y_concave(o));
            x = repmat(x_concave(h),1,length(y));
            y = y';
            x = x';
            
            line = [y,x];
            
            
        elseif slope(h,o) == 0
            x = min(x_concave(h),x_concave(o)):max(x_concave(h),x_concave(o));
            y = repmat(y_concave(h),1,length(x));
            x = x';
            y = y';
            line = [y,x];
            
        else
            x_dist = abs(x_concave(h)-x_concave(o));
            y_dist = abs(y_concave(h)-y_concave(o));
            if x_dist >= y_dist
                x = linspace(min(x_concave(h),x_concave(o)),max(x_concave(h),x_concave(o)),distance([x_concave(h) y_concave(h)], [x_concave(o),y_concave(o)])+1) ;
                x = x';
                for p = 1:length(x(:,1))
                    y(p,1) = slope(h,o)*(x(p,1)-x_concave(h)) + y_concave(h) ;
                end
            elseif y_dist > x_dist
                y = linspace(min(y_concave(h),y_concave(o)),max(y_concave(h),y_concave(o)),distance([x_concave(h) y_concave(h)], [x_concave(o),y_concave(o)])+1) ;
                y = y';
                for p = 1:length(y(:,1))
                    x(p,1) = (y(p,1) - y_concave(h))/slope(h,o) + x_concave(h);
                end
            end
            
            x = floor(x);
            y = floor(y);
            
            line = [y,x];
            
            
        end
        
        y = [];
        x= [];
        
    end
    
end






end