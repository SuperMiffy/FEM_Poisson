function [result] = RHF(x,y)
%right handside function
% input: x
% output: result

% result = pi^2*x.*sin(pi*x)+sin(pi*x)-pi*cos(pi*x);
% result=pi^2/2.*sin(pi*x)-pi./x.*cos(pi*x)-sin(pi*x)./x;
% result=-4;
% result=2*pi^2.*sin(pi.*x).*sin(pi.*y);  
result=2*pi^2.*sin(pi.*x).*sin(pi.*y)+sin(pi.*x).*sin(pi.*y)+pi.*cos(pi.*x).*sin(pi.*y)+pi.*cos(pi.*y).*sin(pi.*x);   

end

