function [result] = x_true_x_der(x,y)
% true solution partial derivative on x-direction

result = pi*cos(pi*x).*sin(pi.*y);
end

