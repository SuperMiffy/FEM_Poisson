function [result] = RHF(x,y,z)
%right handside function
% input: x
% output: result

% result = 3*pi^2.*sin(pi.*x).*sin(pi.*y).*sin(pi.*z);

result=3*pi^2.*sin(pi.*x).*sin(pi.*y).*sin(pi.*z)+pi.*cos(pi.*x).*sin(pi.*y).*sin(pi.*z)+...
    pi.*sin(pi.*x).*cos(pi.*y).*sin(pi.*z)+pi.*sin(pi.*x).*sin(pi.*y).*cos(pi.*z)+...
    sin(pi.*x).*sin(pi.*y).*sin(pi.*z);   
end