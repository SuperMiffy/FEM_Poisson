function [result] = BasisFunc(vertices,x,deg,basis_index,x_der_index,y_der_index,z_der_index)
% generate the basis function on vertices
% input: vertices: start and endpoint of interval
%        x:function variable
%        deg,basis_index,x_der_index,y_der_index,z_der_index
% output: result

a = vertices(:,1);
b = vertices(:,7);

% affine transformation
J11 = b(1)-a(1);
J22 = b(2)-a(2);
J33 = b(3)-a(3);
%J_det = J11*J22*J33;
J_inv = diag([1/J11,1/J22,1/J33]);

x_real = J_inv*([x(:,1)-a(1),x(:,2)-a(2),x(:,3)-a(3)]');%real  node in reference element
x_real = x_real';
if deg == 1
    if x_der_index == 0 && y_der_index == 0 && z_der_index == 0
        result = BasisFuncRefer(x_real,deg,basis_index,0,0,0);
    elseif x_der_index == 1 && y_der_index == 0 && z_der_index == 0
        result = 1/J11*BasisFuncRefer(x_real,deg,basis_index,1,0,0);
    elseif x_der_index == 0 && y_der_index == 1 && z_der_index == 0
        result = 1/J22*BasisFuncRefer(x_real,deg,basis_index,0,1,0);
    elseif x_der_index == 0 && y_der_index == 0 && z_der_index == 1
        result = 1/J33*BasisFuncRefer(x_real,deg,basis_index,0,0,1);
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

function [result] = BasisFuncRefer(x,deg,basis_index,x_der_index,y_der_index,z_der_index)
% generate the basis function on reference vertices
% input: vertices: start and endpoint of interval
%        x:function variable
%        deg,basis_index,x_der_index,y_der_index,z_der_index
% output: result

x_ref = x(:,1);
y_ref = x(:,2);
z_ref = x(:,3);

if deg == 1
    if basis_index == 1
        if x_der_index == 0 && y_der_index == 0 && z_der_index == 0
            result = (1 - x_ref).*(1 - y_ref).*(1 - z_ref);
        elseif x_der_index == 1 && y_der_index == 0 && z_der_index == 0
            result = -1*(1 - y_ref).*(1 - z_ref);
        elseif x_der_index == 0 && y_der_index == 1 && z_der_index == 0
            result = -1*(1 - x_ref).*(1 - z_ref);
        elseif x_der_index == 0 && y_der_index == 0 && z_der_index == 1
            result = -1*(1 - x_ref).*(1 - y_ref);
        end
    elseif basis_index == 2
        if x_der_index == 0 && y_der_index == 0 && z_der_index == 0
            result = x_ref.*(1 - y_ref).*(1 - z_ref);
        elseif x_der_index == 1 && y_der_index == 0 && z_der_index == 0
            result = 1*(1 - y_ref).*(1 - z_ref);
        elseif x_der_index == 0 && y_der_index == 1 && z_der_index == 0
            result = -1*x_ref.*(1 - z_ref);
        elseif x_der_index == 0 && y_der_index == 0 && z_der_index == 1
            result = -1*x_ref.*(1 - y_ref);
        end
    elseif basis_index == 3
        if x_der_index == 0 && y_der_index == 0 && z_der_index == 0
            result = x_ref.*y_ref.*(1 - z_ref);
        elseif x_der_index == 1 && y_der_index == 0 && z_der_index == 0
            result = y_ref.*(1 - z_ref);
        elseif x_der_index == 0 && y_der_index == 1 && z_der_index == 0
            result = x_ref.*(1 - z_ref);
        elseif x_der_index == 0 && y_der_index == 0 && z_der_index == 1
            result = -1*x_ref.*y_ref;
            
        end
    elseif basis_index == 4
        if x_der_index == 0 && y_der_index == 0 && z_der_index == 0
            result = (1 - x_ref).*y_ref.*(1 - z_ref);
        elseif x_der_index == 1 && y_der_index == 0 && z_der_index == 0
            result = -1*y_ref.*(1 - z_ref);
        elseif x_der_index == 0 && y_der_index == 1 && z_der_index == 0
            result = (1 - x_ref).*(1 - z_ref);
        elseif x_der_index == 0 && y_der_index == 0 && z_der_index == 1
            result = -1*(1 - x_ref).*y_ref;
        end
    elseif basis_index == 5
        if x_der_index == 0 && y_der_index == 0 && z_der_index == 0
            result = (1 - x_ref).*(1 - y_ref).*z_ref;
        elseif x_der_index == 1 && y_der_index == 0 && z_der_index == 0
            result = -1*(1 - y_ref).*z_ref;
        elseif x_der_index == 0 && y_der_index == 1 && z_der_index == 0
            result = -1*(1 - x_ref).*z_ref;
        elseif x_der_index == 0 && y_der_index == 0 && z_der_index == 1
            result = (1 - x_ref).*(1 - y_ref);
        end
    elseif basis_index == 6
        if x_der_index == 0 && y_der_index == 0 && z_der_index == 0
            result = x_ref.*(1 - y_ref).*z_ref;
        elseif x_der_index == 1 && y_der_index == 0 && z_der_index == 0
            result = 1*(1 - y_ref).*z_ref;
        elseif x_der_index == 0 && y_der_index == 1 && z_der_index == 0
            result = -1*x_ref.*z_ref;
        elseif x_der_index == 0 && y_der_index == 0 && z_der_index == 1
            result = x_ref.*(1 - y_ref);
        end
    elseif basis_index == 7
        if x_der_index == 0 && y_der_index == 0 && z_der_index == 0
            result = x_ref.*y_ref.*z_ref;
        elseif x_der_index == 1 && y_der_index == 0 && z_der_index == 0
            result = y_ref.*z_ref;
        elseif x_der_index == 0 && y_der_index == 1 && z_der_index == 0
            result = x_ref.*z_ref;
        elseif x_der_index == 0 && y_der_index == 0 && z_der_index == 1
            result = x_ref.*y_ref;
        end
    elseif basis_index == 8
        if x_der_index == 0 && y_der_index == 0 && z_der_index == 0
            result = (1 - x_ref).*y_ref.*z_ref;
        elseif x_der_index == 1 && y_der_index == 0 && z_der_index == 0
            result = -1*y_ref.*z_ref;
        elseif x_der_index == 0 && y_der_index == 1 && z_der_index == 0
            result = (1 - x_ref).*z_ref;
        elseif x_der_index == 0 && y_der_index == 0 && z_der_index == 1
            result = (1 - x_ref).*y_ref;
            
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

