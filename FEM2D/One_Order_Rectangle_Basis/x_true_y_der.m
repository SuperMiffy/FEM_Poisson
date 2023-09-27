function [result] = x_true_y_der(x,y)
% true solution partial derivative on y-direction

result = pi*sin(pi*x).*cos(pi.*y);
end

