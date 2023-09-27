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
Gauss_reference_weights = zeros(nofGauss^2,1);
Gauss_reference_nodes = zeros(nofGauss^2,2);

for i = 1:nofGauss
   Gauss_reference_weights(nofGauss*(i-1)+1:nofGauss*i) = Gauss_1D_reference_weights(i).*...
       Gauss_1D_reference_weights';
   Gauss_reference_nodes(nofGauss*(i-1)+1:nofGauss*i,1) = Gauss_1D_reference_nodes(i);
   Gauss_reference_nodes(nofGauss*(i-1)+1:nofGauss*i,2) = Gauss_1D_reference_nodes';
end
% affine from [-1,1] * [-1,1] into [0,1] * [0,1]
T = [1/2,0;0,1/2];
R = [1/2;1/2];

Gauss_reference_weights = Gauss_reference_weights * abs(det(T));
Gauss_reference_nodes = T*Gauss_reference_nodes'+R;
Gauss_reference_nodes = Gauss_reference_nodes';

% affine into the real element
x1 = vertices(1,1);
y1 = vertices(2,1);
x2 = vertices(1,2);
y2 = vertices(2,2);
x3 = vertices(1,4);
y3 = vertices(2,4);

Jacobi = abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
Gauss_weights = Gauss_reference_weights*Jacobi;
Gauss_nodes(:,1) = (x2-x1).*Gauss_reference_nodes(:,1)+(x3-x1).*Gauss_reference_nodes(:,2)+x1;
Gauss_nodes(:,2) = (y2-y1).*Gauss_reference_nodes(:,1)+(y3-y1).*Gauss_reference_nodes(:,2)+y1;
end

