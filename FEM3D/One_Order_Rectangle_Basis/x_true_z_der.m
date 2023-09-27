function [result] = x_true_z_der(x,y,z)
% true solution partial derivative on z-direction

result = pi*sin(pi*x).*sin(pi.*y).*cos(pi.*z);
end

