function [result] = RHF(x,y)
%right handside function
% input: x
% output: result

% result = 2*pi^2.*sin(pi*x).*sin(pi*y);

result=2*pi^2.*sin(pi.*x).*sin(pi.*y)+sin(pi.*x).*sin(pi.*y)+pi.*cos(pi.*x).*sin(pi.*y)+pi.*cos(pi.*y).*sin(pi.*x);   
end

