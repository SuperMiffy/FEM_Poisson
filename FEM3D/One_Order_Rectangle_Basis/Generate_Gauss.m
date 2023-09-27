function [Gauss_weights,Gauss_nodes] = Generate_Gauss(vertices)
% Gauss integral with five points
% input: vertices
% output: Gauss_weights and Gauss_nodes

% reference element [-1,1] * [-1,1]
Gauss_1D_reference_weights = [0.236926885056189, 0.478628670499366, 0.568888888888889,...
                 0.478628670499366, 0.236926885056189]; 
             
Gauss_1D_reference_nodes = [-0.906179845938664, -0.538469310105683, 0, ...
                0.538469310105683, 0.906179845938664];

nofGauss = length(Gauss_1D_reference_nodes);% number of Gauss integral points
Gauss_reference_weights = zeros(nofGauss^3,1);
Gauss_reference_nodes = zeros(nofGauss^3,3);

for i = 1:nofGauss
    for j = 1:nofGauss
        Gauss_reference_weights(nofGauss^2*(i-1)+nofGauss*(j-1)+1:nofGauss^2*(i-1)+nofGauss*j) = Gauss_1D_reference_weights(i).*...
            Gauss_1D_reference_weights(j).*Gauss_1D_reference_weights';
        Gauss_reference_nodes(nofGauss^2*(i-1)+1:nofGauss^2*i,1) = Gauss_1D_reference_nodes(i);
        Gauss_reference_nodes(nofGauss^2*(i-1)+nofGauss*(j-1)+1:nofGauss^2*(i-1)+nofGauss*j,2) = Gauss_1D_reference_nodes(j);
        Gauss_reference_nodes(nofGauss^2*(i-1)+nofGauss*(j-1)+1:nofGauss^2*(i-1)+nofGauss*j,3) = Gauss_1D_reference_nodes';
    end
end
% affine from [-1,1]^3 into [0,1]^3
T = diag([1/2,1/2,1/2]);
R = [1/2;1/2;1/2];

Gauss_reference_weights = Gauss_reference_weights * abs(det(T));
Gauss_reference_nodes = T*Gauss_reference_nodes'+R;
Gauss_reference_nodes = Gauss_reference_nodes';

% affine into the real element
x1 = vertices(1,1);
y1 = vertices(2,1);
z1 = vertices(3,1);
x7 = vertices(1,7);
y7 = vertices(2,7);
z7 = vertices(3,7);

Jacobi = abs((x7-x1)*(y7-y1)*(z7-z1));
Gauss_weights = Gauss_reference_weights*Jacobi;
Gauss_nodes(:,1) = (x7-x1).*Gauss_reference_nodes(:,1)+x1;
Gauss_nodes(:,2) = (y7-y1).*Gauss_reference_nodes(:,2)+y1;
Gauss_nodes(:,3) = (z7-z1).*Gauss_reference_nodes(:,3)+z1;
end

