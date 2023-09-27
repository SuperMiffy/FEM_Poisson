function [result] = x_true_y_der(x,y,z)
% true solution partial derivative on y-direction

result = pi*sin(pi*x).*cos(pi.*y).*sin(pi.*z);
end

