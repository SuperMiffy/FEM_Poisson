function [result] = RHF(x)
%right handside function
% input: x
% output: result

result = pi^2*x.*sin(pi*x)+sin(pi*x);
% result=pi^2/2.*sin(pi*x)-pi./x.*cos(pi*x)-sin(pi*x)./x;
% result=pi^2/2.*sin(pi*x)-sin(pi*x)./x;
end

