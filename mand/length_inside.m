function [N] = length_inside(c)
%LENGTH_INSIDE Summary of this function goes here
%   Detailed explanation goes here
L = 200;
N = L;
path = zeros(L,1);
path(1) = c;

for k = 2:L
    path(k) = path(k-1).^2 + c;
    if abs(path(k)) > 2
        N = k;
        return
    end
    
end

end

