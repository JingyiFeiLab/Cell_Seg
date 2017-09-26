function new_borders = newBorder(figure);

figure;
field = 'New_Borders';
new_borders = struct(field,[]);
i = 1;
% while 0 < 1
%     [~,~,button] = ginput(1);
%     if isempty(button);
%         break
%     elseif button == 119;
%         h = imfreehand;
%         x = int64(getPosition(h));
%         new_borders(i).New_Borders = x;
%         i = i+1;
%     end
% end
 
h = imfreehand;
x = int64(getPosition(h));
new_borders(i).New_Borders = x;


end
   