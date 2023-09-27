function [Gauss_weights,Gauss_nodes] = Generate_Gauss(vertices)
% Gauss integral with five points
% input: vertices
% output: Gauss_weights and Gauss_nodes

% reference element
Gauss_reference_weights = [64/81*(1-0)/8,100/324*(1-sqrt(3/5))/8,100/324*(1-sqrt(3/5))/8,...
    100/324*(1+sqrt(3/5))/8,100/324*(1+sqrt(3/5))/8,40/81*(1-0)/8,40/81*(1-0)/8,...
    40/81*(1-sqrt(3/5))/8,40/81*(1+sqrt(3/5))/8]';
Gauss_reference_nodes = [(1+0)/2,(1-0)*(1+0)/4;(1+sqrt(3/5))/2,(1-sqrt(3/5))*(1+sqrt(3/5))/4;...
    (1+sqrt(3/5))/2,(1-sqrt(3/5))*(1-sqrt(3/5))/4;(1-sqrt(3/5))/2,(1+sqrt(3/5))*(1+sqrt(3/5))/4;...
    (1-sqrt(3/5))/2,(1+sqrt(3/5))*(1-sqrt(3/5))/4;(1+0)/2,(1-0)*(1+sqrt(3/5))/4;(1+0)/2,(1-0)*(1-sqrt(3/5))/4;...
    (1+sqrt(3/5))/2,(1-sqrt(3/5))*(1+0)/4;(1-sqrt(3/5))/2,(1+sqrt(3/5))*(1+0)/4];

% affine into the real element
x1 = vertices(1,1);
y1 = vertices(2,1);
x2 = vertices(1,2);
y2 = vertices(2,2);
x3 = vertices(1,3);
y3 = vertices(2,3);

Jacobi = abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
Gauss_weights = Gauss_reference_weights*Jacobi;
Gauss_nodes(:,1) = (x2-x1).*Gauss_reference_nodes(:,1)+(x3-x1).*Gauss_reference_nodes(:,2)+x1;
Gauss_nodes(:,2) = (y2-y1).*Gauss_reference_nodes(:,1)+(y3-y1).*Gauss_reference_nodes(:,2)+y1;
end

