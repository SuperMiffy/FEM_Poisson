function [result] = BasisFunc(vertices,x,deg,basis_index,x_der_index,y_der_index)
% generate the basis function on vertices
% input: vertices: start and endpoint of interval
%        x:function variable
%        deg,basis_index,x_der_index,y_der_index
% output: result

a = vertices(:,1);
b = vertices(:,2);
c = vertices(:,3);
d = vertices(:,4);

% affine transformation
J11 = b(1)-a(1);
J12 = d(1)-a(1);
J21 = b(2)-a(2);
J22 = d(2)-a(2);
J_det = det([J11,J12;J21,J22]);
J_inv = [J22,-J12;-J21,J11]/J_det;

x_real = J_inv*([x(:,1)-a(1),x(:,2)-a(2)]');%real  node in reference element
x_real = x_real';
if deg == 1
    if x_der_index == 0 && y_der_index == 0
        result = BasisFuncRefer(x_real,deg,basis_index,0,0);
    elseif x_der_index == 1 && y_der_index == 0
        result = (J22*BasisFuncRefer(x_real,deg,basis_index,1,0)+...
            (-J21)*BasisFuncRefer(x_real,deg,basis_index,0,1))/J_det;
    elseif x_der_index == 0 && y_der_index == 1
        result = ((-J12)*BasisFuncRefer(x_real,deg,basis_index,1,0)+...
            (J11)*BasisFuncRefer(x_real,deg,basis_index,0,1))/J_det;
    elseif x_der_index == 1 && y_der_index == 1
        result = 0;
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

function [result] = BasisFuncRefer(x,deg,basis_index,x_der_index,y_der_index)
% generate the basis function on reference vertices
% input: vertices: start and endpoint of interval
%        x:function variable
%        deg,basis_index,x_der_index,y_der_index
% output: result

x_ref = x(:,1);
y_ref = x(:,2);

if deg == 1
    if basis_index == 1
        if x_der_index == 0 && y_der_index == 0
            result = (1 - x_ref).*(1 - y_ref);
        elseif x_der_index == 1 && y_der_index == 0
            result = y_ref - 1;
        elseif x_der_index == 0 && y_der_index == 1
            result = x_ref - 1;
        elseif x_der_index == 1 && y_der_index == 1
            result = 1;
        end
    elseif basis_index == 2
        if x_der_index == 0 && y_der_index == 0
            result = x_ref.*(1 - y_ref);
        elseif x_der_index == 1 && y_der_index == 0
            result = 1 - y_ref;
        elseif x_der_index == 0 && y_der_index == 1
            result = - x_ref;
        elseif x_der_index == 1 && y_der_index == 1
            result = - 1;
        end
    elseif basis_index == 3
        if x_der_index == 0 && y_der_index == 0
            result = x_ref .* y_ref;
        elseif x_der_index == 1 && y_der_index == 0
            result = y_ref;
        elseif x_der_index == 0 && y_der_index == 1
            result = x_ref;
        elseif x_der_index == 1 && y_der_index == 1
            result = 1;
        end
     elseif basis_index == 4
        if x_der_index == 0 && y_der_index == 0
            result = (1 - x_ref).*y_ref;
        elseif x_der_index == 1 && y_der_index == 0
            result = - y_ref;
        elseif x_der_index == 0 && y_der_index == 1
            result = 1 - x_ref;
        elseif x_der_index == 1 && y_der_index == 1
            result = - 1;
        end
    end
% elseif deg == 2
%     if basis_index == 1
%         if der_index == 0
%             result = 2./(a-b)^2.*(x-(a+b)/2).*(x-b);
%         elseif der_index == 1
%             result = 2./(a-b)^2.*(2.*x-(a+3*b)/2);
%         elseif der_index == 2
%             result = 2*2./(a-b)^2;
%         end
%     elseif basis_index == 2
%         if der_index == 0
%             result = -4./(a-b).^2.*(x-a).*(x-b);
%         elseif der_index == 1
%             result = -4./(a-b).^2.*(2*x-(a+b));
%         elseif der_index == 2
%             result = 2*(-4)./(a-b).^2;
%         end
%     elseif basis_index == 3
%         if der_index == 0
%             result = 2./(a-b)^2.*(x-a).*(x-(a+b)/2);
%         elseif der_index == 1
%             result = 2./(a-b)^2.*(2*x-(3*a+b)/2);
%         elseif der_index == 2
%             result = 2*2./(a-b)^2;
%         end
%     end
end

end

