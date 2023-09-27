function [result] = BasisFunc(vertices,x,basis_index,der_index)
% generate the basis function on vertices
% input: vertices: start and endpoint of interval
%        x:function variable
%        basis_index,der_index
% output: result

a = vertices(1);
b = vertices(2);

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

end

