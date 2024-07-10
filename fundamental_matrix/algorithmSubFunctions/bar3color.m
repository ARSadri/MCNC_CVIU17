function [ ] = bar3color( Mat )
    %Colored bar3
    b = bar3(Mat);
    for k = 1:length(b)
        zdata = get(b(k),'ZData');
        set(b(k),'CData',zdata);
        set(b(k),'FaceColor','interp');
    end
end

