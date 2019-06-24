function [x,y,z] = upd_mand(x0,x1,y0,y1)
x = x0:1/500*(x1-x0):x1;
y = y0:1/500*(y1-y0):y1;

z = zeros(length(x),length(y));

for i = 1:length(x)
    for j = 1:length(y)
        c = x(i) + y(j)*1i;
        z(i,j) = length_inside(c).^-1-1/100;
    end
end

end

