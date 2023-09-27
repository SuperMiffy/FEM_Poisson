function [result] = BasisFunc(vertices,x,deg,basis_index,der_index)
% generate the basis function on vertices
% input: vertices: start and endpoint of interval
%        x:function variable
%        deg,basis_index,der_index
% output: result

a = vertices(1);
b = vertices(2);
if deg == 1
    if basis_index == 1
        if der_index == 0
            result = (b-x)./(b-a);
        elseif der_index == 1
            result = -1/(b-a);
        end
    elseif basis_index == 2
        if der_index == 0
            result = (x-a)./(b-a);
        elseif der_index == 1
            result = 1/(b-a);
        end
    end
elseif deg == 2
    if basis_index == 1
        if der_index == 0
            result = 2./(a-b)^2.*(x-(a+b)/2).*(x-b);
        elseif der_index == 1
            result = 2./(a-b)^2.*(2.*x-(a+3*b)/2);
        elseif der_index == 2
            result = 2*2./(a-b)^2;
        end
    elseif basis_index == 2
        if der_index == 0
            result = -4./(a-b).^2.*(x-a).*(x-b);
        elseif der_index == 1
            result = -4./(a-b).^2.*(2*x-(a+b));
        elseif der_index == 2
            result = 2*(-4)./(a-b).^2;
        end
    elseif basis_index == 3
        if der_index == 0
            result = 2./(a-b)^2.*(x-a).*(x-(a+b)/2);
        elseif der_index == 1
            result = 2./(a-b)^2.*(2*x-(3*a+b)/2);
        elseif der_index == 2
            result = 2*2./(a-b)^2;
        end
    end
end
end

