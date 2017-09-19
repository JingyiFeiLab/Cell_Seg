function BW = smallID(BW)


num = max(BW(:));

for i = 1:num
    %if sum(BWi(:) == i) < 10
    if sum(BW(:) == i) < 10
        %BW(BWi == i) = 0;
        BW(BW == i) = 0;
    end
end

BW;

end