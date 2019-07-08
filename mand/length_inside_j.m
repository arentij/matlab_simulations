function [N] = length_inside_j(z0)
%LENGTH_INSIDE Summary of this function goes here
%   Detailed explanation goes here
L = 100;
N = L;
path = zeros(L,1);
path(1) = z0;
c = -0.4 + 0.6i;
for k = 2:L
    path(k) = path(k-1).^2 + c;
    if abs(path(k)) > 2
        N = k;
        return
    end
    
end

end

